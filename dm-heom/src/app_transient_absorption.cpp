// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include <iostream>
#include <sstream>

#include <boost/format.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <noma/num/meta_stepper.hpp>

#include "heom/command_line.hpp"
#include "heom/common.hpp"
#include "heom/constants.hpp"
#include "heom/dipole_matrix.hpp"
#include "heom/dipole_pseudo_pathway.hpp"
#include "heom/e_field.hpp"
#include "heom/handle_main_exception.hpp"
#include "heom/hierarchy_graph.hpp"
#include "heom/hierarchy_mask.hpp"
#include "heom/hierarchy_norm.hpp"
#include "heom/instance.hpp"
#include "heom/make_file_observer_list.hpp"
#include "heom/population_dynamics_solver.hpp"
#include "heom/ocl_config.hpp"
#include "heom/ode.hpp"
#include "heom/population_dynamics_solver.hpp"
#include "heom/sites_to_states.hpp"
#include "heom/transient_absorption_config.hpp"

namespace bmt = ::heom::bmt;
namespace num = ::heom::num;
namespace ocl = ::heom::ocl;
using heom::int_t;
using heom::real_t;
using heom::real_format;
using heom::complex_t;
using heom::default_delimiter;


int main(int argc, char* argv[])
{
	// start runtime measurement
	bmt::timer app_timer;

	// output compile time configuration
	std::cout << "-------------------- Compile-Time Configuration ------------" << std::endl;
	heom::write_compile_config(std::cout);
	std::cout << "------------------------------------------------------------" << std::endl;

	// command line parsing
	// usage: app_* <ocl_config> <heom_config> [partition_file]
	std::string ocl_config_file = "ocl_config.cfg";
	std::string heom_config_file = "transient_absorption.cfg";

	if (argc >= 2)
		ocl_config_file = argv[1];
	if (argc >= 3)
		heom_config_file = argv[2];

	std::exception_ptr eptr;
	try {
		// command line parsing
		heom::command_line command_line { argc, argv };

		std::cout << "Parsing OpenCL configuration from file: " << command_line.ocl_config_filename() << std::endl;
		heom::ocl_config ocl_config(command_line.ocl_config_filename());
		std::cout << "Parsing HEOM configuration from file: " << command_line.heom_config_filename() << std::endl;
		heom::transient_absorption_config heom_config(command_line.heom_config_filename());

		heom::hierarchy_graph complete_graph(heom_config.baths_number(), heom_config.baths_matsubaras(), heom_config.system_ado_depth());

		using stepper_t = num::meta_stepper; // uses solver_stepper_type from config
		using solver_t = heom::population_dynamics_solver<heom::ode, stepper_t, heom::hierarchy_norm>;

		// prepare output file (so we run into possible file errors before starting the computation)
		auto& observation_filename = heom_config.observations().get().front().get().second; // filename of first specified observation, guaranteed to exist by config checks
		std::ofstream observation_file(observation_filename);
		if (observation_file.fail())
			throw std::runtime_error("Error: could not open output file: " + observation_filename);

		// transient absorption computation over all combinations of dipole matrices and phases
		const size_t num_dipole_matrices = heom_config.dipole_tensor_prefactors().size();
		const size_t num_pulse_phases = heom_config.probe_pulse_phases().size();
		const size_t num_experiments = num_dipole_matrices * num_pulse_phases;

		// solver iterations and steps
		const int_t iterations = heom_config.solver_steps() / heom_config.program_observe_steps();
		const int_t steps_per_iteration = heom_config.program_observe_steps();
		const int_t total_steps = iterations * heom_config.program_observe_steps();

		// output data structure for averaged traces
		const int_t num_observations = iterations + 1;
		std::vector<heom::matrix_trace_observation_view> observations(num_observations);
		// e_field is needed once, and is the same across runs
		std::vector<complex_t> e_field_observations(num_observations);

		// transient absorption pseudo pathway
		const auto& pathway_spec = heom_config.spectra_include_esa() ? heom::dipole_pseudo_pathway_esa_spec : heom::dipole_pseudo_pathway_wgs_spec;
		const auto sts_mode = pathway_spec.sites_to_states_mode();

		// for progress estimate
		const size_t solver_runs = num_dipole_matrices * num_pulse_phases;
		size_t solver_run = 0;

		// iterate over all dipoles and compute population dynamics
		for (size_t tensor_index = 0; tensor_index < num_dipole_matrices; ++tensor_index) {
			// loop over all pulse phases
			for (size_t phase_index = 0; phase_index < num_pulse_phases; ++phase_index)
			{
				const bool first_solver_run = tensor_index == 0 && phase_index == 0; // only observe first solver run

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
				auto get_dipole_matrix = [&](size_t spec_index) {
					return heom::get_dipole_matrix_for_pathway(pathway_spec, spec_index, heom_config.system_sites(), heom_config, tensor_index);
				};
				auto dipole_matrix_minus = get_dipole_matrix(0);
				auto dipole_matrix_plus = get_dipole_matrix(1);
				// TODO: check scaling
				dipole_matrix_minus.scale(heom::get_dipole_tensor_prefactor(pathway_spec, heom_config, tensor_index));

				// TODO: remove this scaling, i.e. fix unit of dipole matrices across all applications, maybe use factor during initialisation of config object or change input file format
				//       scaling dipole strengths in dipole_config_entries does work for this app, but messes with the result of others (tested with linear_absorption)
				dipole_matrix_minus.scale(heom::constants::debye_to_coulomb_meter);
				dipole_matrix_plus.scale(heom::constants::debye_to_coulomb_meter);
//std::cout << dipole_matrix_minus << std::endl << std::endl;
//std::cout << dipole_matrix_plus << std::endl << std::endl;

				// TODO: check if this is really plus and not prescaled... prefactor might be one in config (?) vs. other apps
				// setup observer with pre-scaled dipole matrix
				heom::matrix_trace_observer<solver_t> trace_observer(complete_graph, solver, dipole_matrix_plus.data()); // D+ from left via observer

				// extract hamiltonian as complex_matrix_t for the e_field
				heom::complex_matrix_t hamiltonian(dipole_matrix_plus.rows(), dipole_matrix_plus.cols());
				for(size_t i = 0; i < hamiltonian.rows(); ++i) {
					for(size_t j = 0; j < hamiltonian.cols(); ++j) {
						const size_t real_index = 2 * (i * hamiltonian.rows() + j);
						const size_t imag_index = real_index + 1;
						hamiltonian.at(i, j) = heom::complex_t(heom_instance.hamiltonian()[real_index], heom_instance.hamiltonian()[imag_index]);
					}
				}

				// add e-field pre-evaluation action to the ODE for transient spectra computation
				auto e_field = heom::e_field<decltype(heom_config), heom::ode>(heom_config, phase_index, dipole_matrix_plus, dipole_matrix_minus, hamiltonian);
				solver.ode().clear_pre_evaluate_actions();
				solver.ode().add_pre_evaluate_action(e_field);

				// with create_hierarchy_mask we create mask corresponding to configured filtering strategy, then update solver
				// NOTE: .data() returns pointer to mask_t array
				heom::hierarchy_mask hierarchy_mask = heom::create_hierarchy_mask(heom_config.filtering_strategy(), complete_graph, heom_config.baths_matsubaras(), heom_config.filtering_first_layer());
				solver.update_hierarchy_mask(hierarchy_mask.data());

				std::cout << "-------------------- Hierarchy Mask Counter ----------------" << std::endl;
				heom::write_hierarchy_mask_stats(hierarchy_mask, std::cout);
				std::cout << "------------------------------------------------------------" << std::endl;

				// generate first observation with initial value
				observations[0] += trace_observer.observe_trace(0.0);
				if (first_solver_run)
					e_field_observations[0] = e_field.compute_field(0.0);

				// setup output per propagation output
				heom::config::obervation_type_list per_propagation_observations;

				// build single-propagation observation list from config list
				for (auto& pair : heom_config.observations().get()) {
					auto& obs_type = pair.get().first;
					auto& filename = pair.get().second;

					if (heom::is_single_propagation_observation(obs_type)) {
						std::remove_const<std::remove_reference<decltype(pair)>::type>::type new_pair;
						new_pair.get().first = pair.get().first; // copy observation type

						std::stringstream new_filename;
						size_t dot_pos = filename.find_first_of('.');
						if (dot_pos != std::string::npos) {
							new_filename << filename.substr(0, dot_pos); // copy everything prior to first dot
						} else {
							new_filename << filename; // copy whole filename, in case not dot was found
						}
						new_filename << "_d" << tensor_index << "_p" << phase_index; // add create suffix
						if (dot_pos != std::string::npos) {
							new_filename << filename.substr(dot_pos); // copy file ending
						}

						new_pair.get().second = new_filename.str(); // set new filename
						per_propagation_observations.get().push_back(new_pair);
					}
				}
				heom::observer_list observer_list = heom::make_file_observer_list(per_propagation_observations, complete_graph, solver);
				observer_list.observe(0.0);

				// write file header once
				if (first_solver_run) {
					trace_observer.write_header(observation_file, false);
					observation_file << default_delimiter;
					observation_file << "e-field_real";
					observation_file << default_delimiter;
					observation_file << "e-field_imag";
					observation_file << '\n';
				}

				// propagate system
				for (int_t k = 0; k < iterations; ++k)
				{
					// propagate
					solver.step_forward(steps_per_iteration);

					// after propagation:
					const auto current_step = (k + 1) * steps_per_iteration;
					const real_t current_time = current_step * heom_config.solver_step_size();

					// observe trace of dipole_matrix_minus - rho
					observations[k + 1] += trace_observer.observe_trace(current_time);
					if (first_solver_run)
						e_field_observations[k + 1] = e_field.compute_field(current_time);

					// observe single-propagation observation types
					observer_list.observe(current_time);

					// update status
					std::stringstream status_prefix;
					status_prefix << "Calculation transient absorption ("
					              << "dipole: " << tensor_index + 1 << '/' << num_dipole_matrices
					              << ", "
					              << "phase: " << phase_index + 1 << '/' << num_pulse_phases
					              << "): ";
					heom::write_progress(current_step - 1, total_steps, solver_run, solver_runs, status_prefix.str(), std::cout);
				}

				// write_complex_matrix(result_buffer_top, heom_instance.states(), std::cout);
				std::cout << "-------------------- Solver Runtime Summary ----------------" << std::endl;
				solver.write_runtime_summary(std::cout);
				std::cout << "------------------------------------------------------------" << std::endl;
			}
		}

		// write averaged traces, with added columns for e-field
		auto format = real_format;
		//for (auto& obs : observations) {
		for (int_t i = 0; i < num_observations; ++i) {
			observations[i].avg(num_experiments);
			observation_file << observations[i];
			observation_file << default_delimiter;
			observation_file << format % e_field_observations[i].real();
			observation_file << default_delimiter;
			observation_file << format % e_field_observations[i].imag();
			observation_file << '\n';
		}
		observation_file << std::flush; // tell the OS to write to disk

	} catch (...) {
		eptr = std::current_exception();
	}
	int ret = heom::handle_main_exception(eptr);

	// print application runtime
	std::cout << "main(): application runtime: " << boost::format("%11.2f") % std::chrono::duration_cast<bmt::seconds>(bmt::duration(app_timer.elapsed())).count() << " s" << std::endl;

	return ret;
}
