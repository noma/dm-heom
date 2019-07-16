// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include <exception>
#include <fstream>
#include <iostream>

#include <boost/format.hpp>
#include <noma/num/meta_stepper.hpp>

#include "heom/command_line.hpp"
#include "heom/common.hpp"
#include "heom/dipole_matrix.hpp"
#include "heom/handle_main_exception.hpp"
#include "heom/hierarchy_norm.hpp"
#include "heom/instance.hpp"
#include "heom/make_file_observer_list.hpp"
#include "heom/ocl_config.hpp"
#include "heom/ode.hpp"
#include "heom/population_dynamics_solver.hpp"
#include "heom/two_dimensional_spectra_config.hpp"
#include "heom/sites_to_states.hpp"

namespace bmt = ::heom::bmt;
namespace num = ::heom::num;
namespace ocl = ::heom::ocl;
using heom::int_t;
using heom::real_t;
using heom::complex_t;
using heom::real_format;
using heom::default_delimiter;
using heom::complex_matrix_t;

int main(int argc, char* argv[])
{
	// start runtime measurement
	bmt::timer app_timer;

	// output compile time configuration
	std::cout << "-------------------- Compile-Time Configuration ------------" << std::endl;
	heom::write_compile_config(std::cout);
	std::cout << "------------------------------------------------------------" << std::endl;

	std::exception_ptr eptr;
	try {
		// command line parsing
		heom::command_line command_line { argc, argv };

		std::cout << "Parsing OpenCL configuration from file: " << command_line.ocl_config_filename() << std::endl;
		heom::ocl_config ocl_config(command_line.ocl_config_filename());
		std::cout << "Parsing HEOM configuration from file: " << command_line.heom_config_filename() << std::endl;
		heom::two_dimensional_spectra_config heom_config(command_line.heom_config_filename());

		heom::hierarchy_graph complete_graph(heom_config.baths_number(), heom_config.baths_matsubaras(), heom_config.system_ado_depth());

		using stepper_t = num::meta_stepper; // uses solver_stepper_type from config
		using solver_t = heom::population_dynamics_solver<heom::ode, stepper_t, heom::hierarchy_norm>;

		const int_t steps_per_iteration = heom_config.program_observe_steps();
		const int_t t_1_iterations = heom_config.spectra_steps_t_1() / heom_config.program_observe_steps(); // not including zero
		const int_t t_3_iterations = heom_config.spectra_steps_t_3() / heom_config.program_observe_steps() + 1; // including 0 and max_steps 0 : steps_per_iteration : spectra_t_3_steps

		// prepare output file (so we run into possible file errors before starting the computation)
		auto& observation_filename = heom_config.observations().get().front().get().second; // filename of first specified observation, guaranteed to exist by config checks
		std::ofstream observation_file(observation_filename);
		if (observation_file.fail())
			throw std::runtime_error("Error: could not open output file: " + observation_filename);

		// output data structure for traces
		const int_t num_t_1_observations = t_1_iterations + 1;
		const int_t num_t_3_observations = t_3_iterations;

		// pre-allocate t_3 x t_1 sized 2d vector for accumulating the observations
		//std::vector<std::vector<heom::matrix_trace_observation_view>> observations(num_t_3_observations);
		std::vector<std::vector<complex_t>> observations(num_t_3_observations);
		for (auto& vec : observations) {
			vec.resize(num_t_1_observations, complex_t(0.0, 0.0));
		}

//		std::cout << "Coupling identity: " << heom::state_baths_coupling(heom_config.baths_coupling(), heom::sites_to_states_mode_t::identity) << std::endl;
//		std::cout << "Coupling with ground state: " << heom::state_baths_coupling(heom_config.baths_coupling(), heom::sites_to_states_mode_t::with_ground_state) << std::endl;
//		std::cout << "Coupling excited state absorption: " << heom::state_baths_coupling(heom_config.baths_coupling(), heom::sites_to_states_mode_t::excited_state_absorption) << std::endl;
//
//		std::cout << "Hamiltonian identity: " << heom::state_hamiltonian(heom_config.system_hamiltonian(), heom::sites_to_states_mode_t::identity) << std::endl;
//		std::cout << "Hamiltonian with ground state: " << heom::state_hamiltonian(heom_config.system_hamiltonian(), heom::sites_to_states_mode_t::with_ground_state) << std::endl;
//		std::cout << "Hamiltonian excited state absorption: " << heom::state_hamiltonian(heom_config.system_hamiltonian(), heom::sites_to_states_mode_t::excited_state_absorption) << std::endl;


		const auto& pathways = heom_config.spectra_pathways().get();

		// for progress estimate
		const size_t solver_runs = pathways.size() * heom_config.dipole_tensor_prefactors().size() * t_3_iterations;
		size_t solver_run = 0;

		// iterate over pathways
		for (size_t p = 0; p < pathways.size(); ++p)
		{
			DEBUG_ONLY( std::cout << "Pathway " << p << " is " << pathways.at(p) << std::endl; )
			const auto& pathway_spec = heom::spectra_pathway_to_specification.at(pathways[p]);
			const auto sts_mode = pathway_spec.sites_to_states_mode();

			heom::instance heom_instance(heom_config, sts_mode, complete_graph);
			heom_instance.set_hierarchy_top(); // NOTE: hierarchy top is initialised to zero with first element being complex_t(1.0,0.0) here


			// OpenCL range
			ocl::nd_range range{
				{}, // offset
				{1, static_cast<std::uint64_t>(heom_instance.matrices())}, // global size
				{1, 1} // local size
			};

			// create a solver from configuration and instance
			solver_t solver(ocl_config, range, heom_config, heom_instance);
			std::cout << "-------------------- OpenCL Runtime Configuration ----------" << std::endl;
			solver.write_ocl_runtime_config(std::cout);
			std::cout << "------------------------------------------------------------" << std::endl;

			for (size_t tensor_index = 0; tensor_index < heom_config.dipole_tensor_prefactors().size(); ++tensor_index) {
				for (int_t t_3_index = 0; t_3_index < t_3_iterations; ++t_3_index) {
					const int_t t_3 = t_3_index * steps_per_iteration;

					// re-init instance
					heom_instance.null_hierarchy();
					heom_instance.set_hierarchy_top();
					// re-init solver
					solver.reset();
					++solver_run; // for progress estimate

					// TODO: optimise: multiply top only on host-side (synchronise host/device!) => think through where possible
					auto get_dipole_matrix = [&](size_t spec_index) {
						return heom::get_dipole_matrix_for_pathway(pathway_spec, spec_index, heom_config.system_sites(), heom_config, tensor_index);
					};

					// lambda for generic dipole multiplication with plus/minus, from left/right
					auto dipole_mult = [&](size_t spec_index) {
						const auto dipole_matrix = get_dipole_matrix(spec_index);

						// multiple from left or right on hierarchy according to pathway specification
						if (pathway_spec.pathway()[spec_index].second ==
						    heom::pathway_dipole_side_t::left)
							solver.hierarchy_mmult_left(dipole_matrix.data());
						else // must be right
							solver.hierarchy_mmult_right(dipole_matrix.data());
					};

					// first dipole multiplication from specification
					dipole_mult(0);

					// propagate t_3 steps
					solver.step_forward(t_3);

					// second dipole multiplication from specification
					dipole_mult(1);

					// propagate delay steps
					solver.step_forward(heom_config.spectra_steps_t_delay());

					// third dipole multiplication from specification
					dipole_mult(2);

					// propagate additional t_1 steps and observe

					// setup observer with pre-scaled dipole matrix
					auto dipole_matrix = get_dipole_matrix(3);
					dipole_matrix.scale(heom::get_dipole_tensor_prefactor(pathway_spec, heom_config, tensor_index));
					heom::matrix_trace_observer<solver_t> trace_observer(complete_graph, solver, dipole_matrix.data()); // use last (fourth) polarization for trace, i.e. index 3

					// generate first observation with initial value
					observations[t_3_index][0] += trace_observer.observe_trace(0.0).trace();

					for (int_t t_1_index = 0; t_1_index < t_1_iterations; ++t_1_index) {
						// propagate
						solver.step_forward(steps_per_iteration);

						// after propagation, we are at:
						const auto current_step = (t_1_index + 1) * steps_per_iteration;
						const real_t current_time = current_step * heom_config.solver_step_size();

						// observe
						observations[t_3_index][t_1_index + 1] += trace_observer.observe_trace(current_time).trace();

						// update status
						if(!command_line.no_progress())
							heom::write_progress(current_step - 1, heom_config.spectra_steps_t_1(), solver_run, solver_runs, "Calculating two dimensional spectra: ", std::cout);
					} // t_1
				} // t_3
			} // tensor_index

			// write_complex_matrix(result_buffer_top, heom_instance.states(), std::cout);
			std::cout << "-------------------- Solver Runtime Summary ----------------" << std::endl;
			solver.write_runtime_summary(std::cout);
			std::cout << "------------------------------------------------------------" << std::endl;
		} // pathways


		// TODO: re-enable
		// average traces
//		for (auto& vec : observations)
//			for (auto& elem : vec)
//				elem = { elem.real() / num_dipole_matrices, elem.imag() / num_dipole_matrices };
//				//elem.avg(num_dipole_matrices);

		// write traces
		auto format = real_format;
		observation_file << "t_3" << default_delimiter << "t_1" << default_delimiter << "trace_real" << default_delimiter << "trace_imag" << std::endl;
		for (int_t t_3_index = 0; t_3_index < num_t_3_observations; ++t_3_index) {
			for (int_t t_1_index = 0; t_1_index < num_t_1_observations; ++t_1_index) {
				const auto& current_obs = observations[t_3_index][t_1_index];
				observation_file << format % (t_3_index * steps_per_iteration * heom_config.solver_step_size());
				observation_file << default_delimiter;
				observation_file << format % (t_1_index * steps_per_iteration * heom_config.solver_step_size());
				observation_file << default_delimiter;
				observation_file << format % current_obs.real();
				observation_file << default_delimiter;
				observation_file << format % current_obs.imag();
				observation_file << '\n';
			}
		}
		observation_file << std::flush; // tell the OS to write to disk


		std::cout << "-------------------- Memory Summary ------------------------" << std::endl;
//		std::cout << "local hierarchy buffer size:            "
//		          <<  boost::format("%11.2f MiB\n") % (heom_instance_thermal_state.size_hierarchy_byte() / 1024.0 / 1024.0);
//		std::cout << "local instance size:                    "
//		          <<  boost::format("%11.2f MiB\n") % (heom_instance_thermal_state.allocated_byte() / 1024.0 / 1024.0);
		std::cout << "max. heap via aligned::allocate(..):    "
		          << boost::format("%11.2f MiB") % (heom::memory::aligned::instance().allocated_byte_max() / 1024.0 / 1024.0)
		          << " ("
		          << boost::format("count: %6i") % heom::memory::aligned::instance().allocations_max()
		          << ')' << std::endl;
//		std::cout << "max. OpenCL buffer allocations:         "
//		          <<  boost::format("%11.2f MiB") % (bicgstab.ocl_helper().allocated_byte_max() / 1024.0 / 1024.0)
//		          << " ("
//		          << boost::format("count: %6i") % bicgstab.ocl_helper().allocations_max()
//		          << ')' << std::endl;
//		std::cout << "max. OpenCL buffer allocations:         "
//		          <<  boost::format("%11.2f MiB") % (solver.ocl_helper().allocated_byte_max() / 1024.0 / 1024.0)
//		          << " ("
//		          << boost::format("count: %6i") % solver.ocl_helper().allocations_max()
//		          << ')' << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;

	} catch (...) {
		eptr = std::current_exception();
	}
	int ret = heom::handle_main_exception(eptr);

	// print application runtime
	std::cout << "time in main():\t" << boost::format("%11.2f") % std::chrono::duration_cast<bmt::seconds>(bmt::duration(app_timer.elapsed())).count() << " s" << std::endl;

	return ret;
}


