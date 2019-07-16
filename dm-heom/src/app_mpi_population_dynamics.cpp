// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include <exception>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <noma/num/meta_stepper.hpp>

#include "heom/common.hpp"
#include "heom/communicator.hpp"
#include "heom/distributed_command_line.hpp"
#include "heom/handle_main_exception.hpp"
#include "heom/hierarchy_mask.hpp"
#include "heom/hierarchy_norm.hpp"
#include "heom/hierarchy_partition.hpp"
#include "heom/instance.hpp"
#include "heom/make_file_observer_list.hpp"
#include "heom/ocl_config.hpp"
#include "heom/ode.hpp"
#include "heom/population_dynamics_config.hpp"
#include "heom/population_dynamics_solver.hpp"

namespace bmt = ::heom::bmt;
namespace num = ::heom::num;
namespace ocl = ::heom::ocl;
using num::int_t;
using num::real_t;
using num::real_format;
using heom::rank_t;

int main(int argc, char* argv[])
{
	// start runtime measurement
	bmt::timer app_timer;

	// init communication
	try {
		heom::communicator::init(&argc, &argv);
	} catch (...) {
		heom::handle_main_exception(std::current_exception());
	}
	const rank_t rank = heom::communicator::rank();
	const rank_t ranks = heom::communicator::ranks();

	// output compile time configuration
	if (rank == 0) {
		std::cout << "-------------------- Compile-Time Configuration ------------" << std::endl;
		heom::write_compile_config(std::cout);
		std::cout << "------------------------------------------------------------" << std::endl;
	}

	std::exception_ptr eptr;
	try {
		// command line parsing
		heom::distributed_command_line command_line { argc, argv };

		std::cout << "main(..), rank " << rank << ": parsing OpenCL configuration from file: " << command_line.ocl_config_filename() << std::endl;
		heom::ocl_config ocl_config(command_line.ocl_config_filename());
		std::cout << "main(..), rank " << rank << ": parsing HEOM configuration from file: " << command_line.heom_config_filename() << std::endl;
		heom::population_dynamics_config heom_config(command_line.heom_config_filename());

		// create hierarchy graph and partitioning
		heom::hierarchy_graph complete_graph(heom_config.baths_number(), heom_config.baths_matsubaras(), heom_config.system_ado_depth());

		heom::partition_mapping partitioning;
		if (command_line.partition_mapping_filename().empty()) {
			std::cout << "main(..), rank " << rank << ": calculating partitioning..." << std::endl;
			// partitioning = create_simple_partitioning(complete_graph, rank, ranks);
			// partitioning = create_round_robin_partitioning(complete_graph, rank, ranks);
			partitioning = create_layer_partitioning(complete_graph, rank, ranks);
		} else {
			std::cout << "main(..), rank " << rank << ": creating partitioning from file '" << command_line.partition_mapping_filename() << "'..." << std::endl;
			partitioning = create_file_partitioning(complete_graph, rank, ranks, command_line.partition_mapping_filename());
		}

		const rank_t master_rank = partitioning.at(0); // NOTE: whoever owns the hierarchy top is the master rank and generates the observation

		heom::hierarchy_partition graph_partition(complete_graph, partitioning, rank, ranks);

		DEBUG_ONLY( std::cout << "main(..), rank " << rank << ": complete_graph.first_non_plus_id(), uid = " << complete_graph.first_non_plus_uid() << ", graph_partition.first_non_plus_id(), lid = " << graph_partition.first_non_plus_uid() << ", uid = " << graph_partition.lid_to_uid().at(graph_partition.first_non_plus_uid()) << std::endl; )

		heom::instance heom_instance(heom_config, heom::sites_to_states_mode_t::identity, graph_partition);

		heom::communicator comm(graph_partition, heom_instance);

		if (rank == master_rank) {
			heom_instance.set_hierarchy_top(heom_config.population_dynamics_rho_init().data());
		}

		using stepper_t = num::meta_stepper; // uses solver_stepper_type from config
		using solver_t = heom::population_dynamics_solver<heom::ode, stepper_t, heom::hierarchy_norm>;

		// OpenCL range
		ocl::nd_range range {
			{}, // offset
			{1, static_cast<std::uint64_t>(graph_partition.compute_nodes())}, // static_cast<std::uint64_t>(heom_instance.matrices())}, // global size
			{1, 1} // local size
		};

		// create solver
		solver_t solver(ocl_config, range, heom_config, heom_instance);
		if (rank == master_rank) {
			std::cout << "-------------------- OpenCL Runtime Configuration ----------" << std::endl;
			solver.write_ocl_runtime_config(std::cout);
			std::cout << "------------------------------------------------------------" << std::endl;
		}

		// with create_hierarchy_mask we create mask corresponding to configured filtering strategy, then update solver
		// NOTE: .data() returns pointer to mask_t array
		heom::hierarchy_mask hierarchy_mask = heom::create_hierarchy_mask(heom_config.filtering_strategy(), complete_graph, heom_config.baths_matsubaras(), heom_config.filtering_first_layer());
		solver.update_hierarchy_mask(hierarchy_mask.data());

		std::cout << "-------------------- Hierarchy Mask Counter ----------------" << std::endl;
		heom::write_hierarchy_mask_stats(hierarchy_mask, std::cout);
		std::cout << "------------------------------------------------------------" << std::endl;

		// register neighbour exchange before each ODE evaluation
		using ne_action_t = heom::neighbour_exchange_action<heom::ode>;
		ne_action_t ne_action(comm, heom_instance);
		solver.ode().add_pre_evaluate_action(std::ref(ne_action));

		heom::observer_list observer_list;

		// boost::format object to re-use for output
		auto format = real_format;

		if (rank == master_rank) {
			// setup output
			observer_list = heom::make_file_observer_list(heom_config, complete_graph, solver);
			// generate first observation with initial value
			observer_list.observe(0.0);
		}

		// run the solver, print output every n steps
		int_t iterations = heom_config.solver_steps() / heom_config.program_observe_steps();
		int_t steps_per_iteration = heom_config.program_observe_steps();
		int_t total_steps = iterations * heom_config.program_observe_steps();

		for (int_t i = 0; i < iterations; ++i) {
			// propagate
			solver.step_forward(steps_per_iteration);

			// after propagation, we are at:
			const auto current_step = (i + 1) * steps_per_iteration;
			const real_t current_time = current_step * heom_config.solver_step_size();

			if (rank == master_rank) {
				// observe
				observer_list.observe(current_time);

				// update status
				heom::write_progress(current_step - 1, total_steps, "Calculation population dynamics: ", std::cout);
			}
		}

		// create sub-directory for profiling data inside current working directory
		boost::filesystem::path profiling_dir { "profile" };
		if (rank == master_rank) { // assuming shared filesystem
			if (!boost::filesystem::create_directory(profiling_dir)) {
				std::cout << "Warning: cannot create profiling output directory: " << profiling_dir << ", existing data might be overwritten." << std::endl;
			}
		}
		comm.barrier(); // make sure no rank is ahead of directory creation on master rank

		// globally collect profiling data for solver across all ranks on the master rank
		heom::statistics_vector solver_stats { solver.stats() }; // NOTE: implied order: step is expected come first
		std::vector<std::vector<heom::runtime_summary_record>> solver_records; // outer vector has all different statistics, inner vector all ranks
		for (const auto& current_pair : solver_stats)
			solver_records.push_back(comm.gather(master_rank, heom::runtime_summary_record(rank, current_pair.first))); // first is statistics object

		// globally collect profiling data for neighbour exchange action across all ranks on the master rank
		heom::statistics_vector ne_action_stats { ne_action.stats() };
		std::vector<std::vector<heom::runtime_summary_record>> ne_action_records; // outer vector has all different statistics, inner vector all ranks
		for (const auto& current_pair : ne_action_stats)
			ne_action_records.push_back(comm.gather(master_rank, heom::runtime_summary_record(rank, current_pair.first))); // first is statistics object

		// globally collect profiling data memory
		std::vector<std::string> memory_names { "local instance size",
		                                        "  hierarchy buffer size",
		                                        "    computed hierarchy nodes",
		                                        "    halo hierarchy nodes",
		                                        "max. heap via aligned::allocate(..)",
		                                        "max. OpenCL buffer allocations" };
		// NOTE: order must match names
		std::vector<std::vector<heom::memory_summary_record>> memory_records; // outer vector has all different statistics, inner vector all ranks
		memory_records.push_back(comm.gather(master_rank, heom::memory_summary_record(rank, 1, heom_instance.allocated_byte())));
		memory_records.push_back(comm.gather(master_rank, heom::memory_summary_record(rank, 1, heom_instance.size_hierarchy_byte()))); // first is statistics object
		memory_records.push_back(comm.gather(master_rank, heom::memory_summary_record(rank, 1, heom_instance.size_hierarchy_top_byte() * graph_partition.compute_nodes())));
		memory_records.push_back(comm.gather(master_rank, heom::memory_summary_record(rank, 1, heom_instance.size_hierarchy_top_byte() * graph_partition.halo_nodes())));
		memory_records.push_back(comm.gather(master_rank, heom::memory_summary_record(rank, heom::memory::aligned::instance().allocations_max(), heom::memory::aligned::instance().allocated_byte_max())));
		memory_records.push_back(comm.gather(master_rank, heom::memory_summary_record(rank, solver.ocl_helper().allocations_max(), solver.ocl_helper().allocated_byte_max())));

		if (rank == master_rank) {
			// find rank with best and worst solver_step
			rank_t best_rank = 0;
			rank_t worst_rank = 0;
			auto& solver_step_records = solver_records[0]; // NOTE: implied order
			for (rank_t r = 0; r < ranks; ++r) {
				best_rank = solver_step_records[r].sum() < solver_records[0][best_rank].sum() ? r : best_rank; // find rank with minimum sum
				worst_rank = solver_step_records[r].sum() > solver_records[0][worst_rank].sum() ? r : worst_rank; // find rank with minimum sum
			}

			// create data structure for and write the solver summary table
			std::vector<heom::runtime_summary_table_entry> solver_table_entries;
			for (size_t i = 0; i < solver_stats.size(); ++i)
				solver_table_entries.emplace_back(solver_stats[i].second, solver_records[i][best_rank], solver_records[i][worst_rank]); // second is the name of the record
			heom::write_runtime_summary_table("Solver Runtime Summary", solver_table_entries, std::cout);

			// create data structure for and write the neighbour exchange summary table
			std::vector<heom::runtime_summary_table_entry> ne_action_table_entries;
			for (size_t i = 0; i < ne_action_stats.size(); ++i)
				ne_action_table_entries.emplace_back(ne_action_stats[i].second, ne_action_records[i][best_rank], ne_action_records[i][worst_rank]); // second is the name of the record
			heom::write_runtime_summary_table("Neighbour Exchange Runtime Summary", ne_action_table_entries, std::cout);

			// create data structure for and write the memory summary table
			std::vector<heom::memory_summary_table_entry> memory_table_entries;
			for (size_t i = 0; i < memory_names.size(); ++i)
				memory_table_entries.emplace_back(memory_names[i], memory_records[i][best_rank], memory_records[i][worst_rank]);
			heom::write_memory_summary_table("Memory Allocation Summary", memory_table_entries, std::cout);

			// write communication matrix
			heom::write_to_file((profiling_dir / "neighbour_exchange_matrix.dat").string(), std::bind(&heom::communicator::write_neighbour_exchange_table, &comm, std::placeholders::_1)); // comm.write_neighbour_exchange_table(neighbour_exchange_matrix_file)
			heom::write_to_file((profiling_dir / "neighbour_exchange_matrix_raw.dat").string(), std::bind(&heom::communicator::write_neighbour_exchange_table_raw, &comm, std::placeholders::_1)); // comm.write_neighbour_exchange_table_raw(neighbour_exchange_matrix_file)
		}

		// output detailed data for every rank separately
		const std::string runtime_stats_filename { "solver_statistics_rank_" + boost::lexical_cast<std::string>(rank) + ".dat" };
		heom::write_to_file((profiling_dir / runtime_stats_filename).string(), std::bind(&solver_t::write_runtime_stats, &solver, std::placeholders::_1)); // solver.write_runtime_stats(runtime_stats_file));

		const std::string ne_action_stats_filename { "neighbour_exchange_statistics_rank_" + boost::lexical_cast<std::string>(rank) + ".dat" };
		heom::write_to_file((profiling_dir / ne_action_stats_filename).string(), std::bind(&ne_action_t::write_runtime_stats, &ne_action, std::placeholders::_1)); // solver.write_runtime_stats(runtime_stats_file));
		
		// copy /proc/curproc/status file for each rank into profiling dir, if exists
		pid_t pid = getpid();
		boost::filesystem::path status_path { "/proc/" + boost::lexical_cast<std::string>(pid) + "/status" };
		boost::filesystem::path status_copy_path { profiling_dir / ("proc_status_rank_" + boost::lexical_cast<std::string>(rank)) };
		if (boost::filesystem::exists(status_path)) {
			if (boost::filesystem::exists(status_copy_path))
				boost::filesystem::remove(status_copy_path); // delete if exists, copy_file line above with overwrite_if_exists has permission issues
			//bost::filesystem::copy_file(status_path, status_copy_path, boost::filesystem::copy_option::overwrite_if_exists);
			//boost::filesystem::copy(status_path, status_copy_path);
			//boost::filesystem::permissions(status_copy_path, boost::filesystem::perms::add_perms | boost::filesystem::perms::owner_write); // make sure file has write rights, as source is read-only
			std::ifstream status_file(status_path.string());
			std::ofstream status_copy_file(status_copy_path.string());
			status_copy_file << status_file.rdbuf();
		}
	} catch (...) {
		eptr = std::current_exception();
	}

	// deinit communication
	try {
		heom::communicator::deinit();
	} catch (...) {
		heom::handle_main_exception(std::current_exception());
	}
	int ret = heom::handle_main_exception(eptr);

	// print application runtime
	if (rank == 0) // NOTE: master_rank out of scope already
		std::cout << "main(): rank " << rank << ": application runtime: " << boost::format("%11.2f") % std::chrono::duration_cast<bmt::seconds>(bmt::duration(app_timer.elapsed())).count() << " s" << std::endl;

	return ret;
}
