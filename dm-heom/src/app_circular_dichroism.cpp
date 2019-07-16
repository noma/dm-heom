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
#include "heom/dipole_pseudo_pathway.hpp"
#include "heom/handle_main_exception.hpp"
#include "heom/hierarchy_mask.hpp"
#include "heom/hierarchy_norm.hpp"
#include "heom/instance.hpp"
#include "heom/make_file_observer_list.hpp"
#include "heom/ocl_config.hpp"
#include "heom/ode.hpp"
#include "heom/circular_dichroism_config.hpp"
#include "heom/population_dynamics_solver.hpp"
#include "heom/sites_to_states.hpp"

namespace bmt = ::heom::bmt;
namespace num = ::heom::num;
namespace ocl = ::heom::ocl;
using num::int_t;
using num::real_t;

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
		heom::circular_dichroism_config heom_config(command_line.heom_config_filename());

		heom::hierarchy_graph complete_graph(heom_config.baths_number(), heom_config.baths_matsubaras(), heom_config.system_ado_depth());

		using stepper_t = num::meta_stepper; // uses solver_stepper_type from config
		using solver_t = heom::population_dynamics_solver<heom::ode, stepper_t, heom::hierarchy_norm>;

		// prepare output file (so we run into possible file errors before starting the computation)
		auto& observation_filename = heom_config.observations().get().front().get().second; // filename of first specified observation, guaranteed to exist by config checks
		std::ofstream observation_file(observation_filename);
		if (observation_file.fail())
			throw std::runtime_error("Error: could not open output file: " + observation_filename);

		int_t iterations = heom_config.solver_steps() / heom_config.program_observe_steps();
		int_t steps_per_iteration = heom_config.program_observe_steps();
		int_t total_steps = iterations * heom_config.program_observe_steps();

		//output data structure for traces
		const size_t num_dipole_matrices = heom_config.dipole_tensor_prefactors().size();
		const int_t num_observations = iterations + 1;
		std::vector<heom::matrix_trace_observation_view> observations(num_observations);

		// circular dichroism pseudo pathway with ground state
		const auto& pathway_spec = heom::dipole_pseudo_pathway_wgs_spec;
		const auto sts_mode = pathway_spec.sites_to_states_mode();

		// for progress estimate
		const size_t solver_runs = num_dipole_matrices;
		size_t solver_run = 0;

		// iterate over all dipoles and compute population dynamics
		for (size_t tensor_index = 0; tensor_index < heom_config.dipole_tensor_prefactors().size(); ++tensor_index) {
			heom::instance heom_instance(heom_config, sts_mode, complete_graph);
			heom_instance.set_hierarchy_top(); // NOTE: hierarchy top is initialised to zero with first element being complex_t(1.0,0.0) here

			// OpenCL range
			ocl::nd_range range {
				{}, // offset
				{1, static_cast<std::uint64_t>(heom_instance.matrices())}, // global size
				{1, 1} // local size
			};

			// create a solver from configuration and instance
			solver_t solver(ocl_config, range, heom_config, heom_instance);
			std::cout << "-------------------- OpenCL Runtime Configuration ----------" << std::endl;
			solver.write_ocl_runtime_config(std::cout);
			std::cout << "------------------------------------------------------------" << std::endl;
			++solver_run; // for progress estimate

			// compute dipole matrices
			auto dipole_matrix_minus = heom::get_dipole_matrix_for_pathway_circular(pathway_spec, 0, heom_config.system_sites(), heom_config, tensor_index);
			auto dipole_matrix_plus  = heom::get_dipole_matrix_for_pathway(pathway_spec, 1, heom_config.system_sites(), heom_config, tensor_index);

			// pre-scale D- for observer
			dipole_matrix_minus.scale(heom::get_dipole_tensor_prefactor(pathway_spec, heom_config, tensor_index));

			// multiply D+ from left
			solver.hierarchy_mmult_left(dipole_matrix_plus.data());

			// setup observer with pre-scaled D-
			heom::matrix_trace_observer<solver_t> trace_observer(complete_graph, solver, dipole_matrix_minus.data()); // D- from left via observer

			// generate first observation with initial value
			observations[0] += trace_observer.observe_trace(0.0);

			// write header once
			if (tensor_index == 0)
				trace_observer.write_header(observation_file, true);

			for (int_t i = 0; i < iterations; ++i) {
				// propagate
				solver.step_forward(steps_per_iteration);

				// after propagation, we are at:
				const auto current_step = (i + 1) * steps_per_iteration;
				const real_t current_time = current_step * heom_config.solver_step_size();

				// observe
				observations[i + 1] += trace_observer.observe_trace(current_time);

				// update status
				heom::write_progress(current_step - 1, total_steps, solver_run, solver_runs, "Calculation circular dichroism: ", std::cout);
			}

			// write_complex_matrix(result_buffer_top, heom_instance.states(), std::cout);
			std::cout << "-------------------- Solver Runtime Summary ----------------" << std::endl;
			solver.write_runtime_summary(std::cout);
			std::cout << "------------------------------------------------------------" << std::endl;

		} // for tensor_index

		// write averaged traces
		for (auto& obs : observations) {
			obs.avg(num_dipole_matrices);
			observation_file << obs << '\n';
		}
		observation_file << std::flush; // tell the OS to write to disk

		std::cout   << "-------------------- Memory Summary ------------------------" << std::endl;
		std::cout   << "max. heap via aligned::allocate(..):    "
		            << boost::format("%11.2f MiB") % (heom::memory::aligned::instance().allocated_byte_max() / 1024.0 / 1024.0)
		            << " ("
		            << boost::format("count: %6i") % heom::memory::aligned::instance().allocations_max()
		            << ')' << std::endl;
//		std::cout   << "max. OpenCL buffer allocations:         "
//		            << boost::format("%11.2f MiB") % (solver.ocl_helper().allocated_byte_max() / 1024.0 / 1024.0)
//		            << " ("
//		            << boost::format("count: %6i") % solver.ocl_helper().allocations_max()
//		            << ')' << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;

	} catch (...) {
		eptr = std::current_exception();
	}
	int ret = heom::handle_main_exception(eptr);

	// print application runtime
	std::cout << "time in main():\t" << boost::format("%11.2f") % std::chrono::duration_cast<bmt::seconds>(bmt::duration(app_timer.elapsed())).count() << " s" << std::endl;

	return ret;
}


