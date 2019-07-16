// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include <algorithm>
#include <exception>
#include <iostream>

#include <boost/format.hpp>
#include <noma/num/meta_stepper.hpp>

#include "heom/command_line.hpp"
#include "heom/common.hpp"
#include "heom/handle_main_exception.hpp"
#include "heom/hierarchy_graph.hpp"
#include "heom/hierarchy_mask.hpp"
#include "heom/hierarchy_norm.hpp"
#include "heom/instance.hpp"
#include "heom/make_file_observer_list.hpp"
#include "heom/ocl_config.hpp"
#include "heom/ode.hpp"
#include "heom/population_dynamics_solver.hpp"
#include "heom/rate_kernel_config.hpp"

namespace bmt = ::heom::bmt;
namespace num = ::heom::num;
namespace ocl = ::heom::ocl;
using num::int_t;
using num::real_t;
using num::complex_t;
using heom::complex_matrix_t;

int main(int argc, char* argv[])
{
	// start runtime measurement
	bmt::timer app_timer;

	// output compile-time configuration
	std::cout << "-------------------- Compile-Time Configuration ------------" << std::endl;
	heom::write_compile_config(std::cout);
	std::cout << "------------------------------------------------------------" << std::endl;

	std::exception_ptr eptr;
	try {
		// command line parsing
		heom::command_line command_line { argc, argv };

		// parse config files
		std::cout << "Parsing OpenCL configuration from file: " << command_line.ocl_config_filename() << std::endl;
		heom::ocl_config ocl_config(command_line.ocl_config_filename());
		std::cout << "Parsing HEOM configuration from file: " << command_line.heom_config_filename() << std::endl;
		heom::rate_kernel_config heom_config(command_line.heom_config_filename());

		heom::hierarchy_graph complete_graph(heom_config.baths_number(), heom_config.baths_matsubaras(), heom_config.system_ado_depth());

		heom::instance heom_instance(heom_config, heom::sites_to_states_mode_t::identity, complete_graph);
		heom_instance.set_hierarchy_top(); // NOTE: default hierarchy_top, i.d. upper left = (1.0, 0.0)

		using stepper_t = num::meta_stepper; // uses solver_stepper_type from config
		using solver_t = heom::population_dynamics_solver<heom::ode, stepper_t, heom::hierarchy_norm>;

		// create a solver from configuration and instance
		const ocl::nd_range ocl_nd_range {
				{}, // offset
		        {1, static_cast<std::uint64_t>(heom_instance.matrices())}, // global size
		        {1, 1} // local size
			};

		solver_t solver(ocl_config, ocl_nd_range, heom_config, heom_instance);
		std::cout << "-------------------- OpenCL Runtime Configuration ----------" << std::endl;
		solver.write_ocl_runtime_config(std::cout);
		std::cout << "------------------------------------------------------------" << std::endl;

		// with create_hierarchy_mask we create mask corresponding to configured filtering strategy, then update solver
		// NOTE: .data() returns pointer to mask_t array
		heom::hierarchy_mask hierarchy_mask = heom::create_hierarchy_mask(heom_config.filtering_strategy(), complete_graph, heom_config.baths_matsubaras(), heom_config.filtering_first_layer());
		solver.update_hierarchy_mask(hierarchy_mask.data());

		std::cout << "-------------------- Hierarchy Mask Counter ----------------" << std::endl;
		heom::write_hierarchy_mask_stats(hierarchy_mask, std::cout);
		std::cout << "------------------------------------------------------------" << std::endl;

		// hierarchy top buffer, needed throughout the computation
		std::unique_ptr<complex_t[]> hierarchy_top_buffer(new complex_t[heom_instance.size_hierarchy_top()]);

		// lambda to evaluate the HEOM ODE outside of the solver's propagation, i.e. non-intrusively (only back_buffer will be overwritten)
		// TODO: think: heom_kernel needs limited depth, probably at most 2 layers
		auto evaluate_heom_ode = [&]() {
			solver.ode().solve(solver.time(), 1.0, solver.front_buffer(), solver.back_buffer()); // NOTE: back_buffer() can always be safely overwritten
		};

		// lambda to scale the whole hierarchy in back_buffer with i, store result to front_buffer
		auto scale_hierarchy_with_i = [&]() {
			solver.hierarchy_scale(solver.back_buffer(), complex_t(0.0,1.0), solver.front_buffer());
		};

		// set hierarchy_top diagonal to zero in front_buffer
		auto zero_hierarchy_top_diag = [&]() {
			// read hierarchy_top from back_buffer
			solver.hierarchy_top(hierarchy_top_buffer.get());
			// modify diagonal
			for (int_t i = 0; i < heom_config.system_sites(); ++i)
				hierarchy_top_buffer[i * heom_config.system_sites() + i] = complex_t { 0.0, 0.0 };
			// update_hierarchy_top
			solver.update_hierarchy_top(hierarchy_top_buffer.get());
		};

		// lambda to read the diagonal of the top level hierarchy
		auto read_top_level_diag = [&](int_t site, complex_matrix_t& matrix_buffer) {
			// read top matrix
			solver.read_buffer(solver.back_buffer(), hierarchy_top_buffer.get(), heom_instance.size_hierarchy_top());
			// copy diagonal to site-th row
			for (int_t i = 0; i < heom_config.system_sites(); ++i)
				matrix_buffer.at(site, i) = hierarchy_top_buffer[i * heom_config.system_sites() + i];
		};


		int_t iterations = heom_config.solver_steps() / heom_config.program_observe_steps();
		int_t steps_per_iteration = heom_config.program_observe_steps();
		int_t total_steps = iterations * heom_config.program_observe_steps();

		// set up output file
		auto& observation_filename = heom_config.observations().get().front().get().second; // filename of first specified observation, guaranteed to exist by config checks
		std::ofstream observation_file(observation_filename);
		if (observation_file.fail())
			throw std::runtime_error("Error: could not open output file: " + observation_filename);

		// set up observation data structures
		std::vector<real_t> observation_times;
		std::vector<complex_matrix_t> observations; // for collecting all observations first
		for (int_t i = 0; i <= iterations; ++i) // add iterations + 1 matrices
			observations.emplace_back(heom_config.system_sites(), heom_config.system_sites()); // NOTE: null-initialised

		// for progress estimate
		const size_t solver_runs = heom_config.system_sites();
		size_t solver_run = 1;

		// iterate over the number of sites
		for (int_t site = 0; site < heom_config.system_sites(); ++site) {
			std::cout << "Processing site " << (site + 1) << "/" << heom_config.system_sites() << std::endl;

			// re-use solver and instance instead of creating a new one for every iteration
			if (site > 0) {
				// re-init instance
				heom_instance.null_hierarchy();

				// create and hierarchy_top
				std::fill(hierarchy_top_buffer.get(), hierarchy_top_buffer.get() + heom_instance.size_hierarchy_top(), complex_t {0.0, 0.0});
				hierarchy_top_buffer[site * heom_config.system_sites() + site] = complex_t { 1.0, 0.0 };
				heom_instance.set_hierarchy_top(hierarchy_top_buffer.get());

				// re-init solver
				solver.reset();
				++solver_run; // for progress estimate
			}


			// run heom_kernel (front_buffer => back_buffer)
			evaluate_heom_ode();
			// multiply with i (imag) (back_buffer => front_buffer)
			scale_hierarchy_with_i();
			// zero diagonal (front_buffer => front_buffer)
			zero_hierarchy_top_diag();
			// run heom_kernel (front_buffer => back_buffer)
			evaluate_heom_ode();
			// observe back_buffer
			read_top_level_diag(site, observations[0]);
			observation_times.push_back(0.0);

			// run the solver, print output every n steps
			for (int_t i = 0; i < iterations; ++i) {
				// propagate
				solver.step_forward(steps_per_iteration); // NOTE: switches front_buffer and back_buffer, result is new front_buffer
				// zero diagonal (front_buffer => front_buffer)
				zero_hierarchy_top_diag();
				// run heom_kernel (front_buffer => back_buffer)
				evaluate_heom_ode();
				// observe back_buffer
				read_top_level_diag(site, observations[i + 1]);

				// after propagation, we are at:
				const auto current_step = (i + 1) * steps_per_iteration;
				const real_t current_time = current_step * heom_config.solver_step_size();

				observation_times.push_back(current_time);

				// update status
				heom::write_progress(current_step - 1, total_steps, solver_run, solver_runs, "Calculation rate kernel: ", std::cout);
			}

		} // for system_sites

		// write observations header
		observation_file << "time";
		for (int_t i = 0; i < heom_config.system_sites(); ++i)
			for (int_t j = 0; j < heom_config.system_sites(); ++j)
				for (int_t k = 0; k < 2; ++k) // complex
				{
					observation_file << heom::default_delimiter;
					observation_file << "elem_" << i << "_" << j << (k % 2 == 0 ? "_real" : "_imag");
				}
		observation_file << std::endl;

		// write observations data
		auto format_0 = heom::real_format;
		auto format_1 = heom::real_format;
		for (int_t it = 0; it <= iterations; ++it) {
			observation_file << (format_0 % observation_times[it]) << heom::default_delimiter;
			for (int_t i = 0; i < heom_config.system_sites(); ++i)
				for (int_t j = 0; j < heom_config.system_sites(); ++j) {
					observation_file << (format_0 % observations[it].at(i, j).real()) << heom::default_delimiter << (format_1 % observations[it].at(i, j).imag());
					if ((i * heom_config.system_sites() + j) < (heom_config.system_sites() * heom_config.system_sites()) - 1)
						observation_file << heom::default_delimiter;
				}
			observation_file << '\n';
		}
		observation_file << std::flush; // tell the OS to write to disk

		// output benchmark results
		std::cout << "-------------------- Solver Runtime Summary ----------------" << std::endl;
		solver.write_runtime_summary(std::cout);
		std::cout << "------------------------------------------------------------" << std::endl;
		std::cout << "-------------------- Memory Summary ------------------------" << std::endl;
		std::cout << "local hierarchy buffer size:            "
		          << boost::format("%11.2f MiB\n") % (heom_instance.size_hierarchy_byte() / 1024.0 / 1024.0);
		std::cout << "local instance size:                    "
		          << boost::format("%11.2f MiB\n") % (heom_instance.allocated_byte() / 1024.0 / 1024.0);
		std::cout << "max. heap via aligned::allocate(..):    "
		          << boost::format("%11.2f MiB") %
		             (heom::memory::aligned::instance().allocated_byte_max() / 1024.0 / 1024.0)
		          << " ("
		          << boost::format("count: %6i") % heom::memory::aligned::instance().allocations_max()
		          << ')' << std::endl;
		std::cout << "max. OpenCL buffer allocations:         "
		          << boost::format("%11.2f MiB") % (solver.ocl_helper().allocated_byte_max() / 1024.0 / 1024.0)
		          << " ("
		          << boost::format("count: %6i") % solver.ocl_helper().allocations_max()
		          << ')' << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;

	} catch (...) {
		eptr = std::current_exception();
	}
	int ret = heom::handle_main_exception(eptr);

	// print application runtime
	std::cout << "main(): application runtime: " << boost::format("%11.2f") % std::chrono::duration_cast<bmt::seconds>(bmt::duration(app_timer.elapsed())).count() << " s" << std::endl;

	return ret;
}
