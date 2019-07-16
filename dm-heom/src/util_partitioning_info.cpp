// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include <exception>
#include <iostream>

#include <boost/format.hpp>

#include "heom/common.hpp"
#include "heom/handle_main_exception.hpp"
#include "heom/hierarchy_partition.hpp"
#include "heom/population_dynamics_config.hpp"

using heom::rank_t;

int main(int argc, char* argv[])
{
	// command line arguments, positional
	// usage:
	//     util_graph_to_file <heom_config_file> <ranks> [partition_file]
	std::string heom_config_file = "population_dynamics.cfg";
	rank_t ranks = 0;
	std::string partition_file = ""; //heom_config_file + ".graph.part." + boost::lexical_cast<std::string>(ranks);

	//if (argc >= 2)
	heom_config_file = argv[1];
	//if (argc >= 3)
	ranks = boost::lexical_cast<rank_t>(argv[2]);
	if (argc >= 4)
		partition_file = argv[3];

	std::exception_ptr eptr;
	try {
		heom::population_dynamics_config heom_config(heom_config_file);

		// create hierarchy graph
		std::cout << "main(..): creating hierarchy_graph..." << std::endl;
		heom::hierarchy_graph complete_graph(heom_config.baths_number(), heom_config.baths_matsubaras(), heom_config.system_ado_depth());

		// create partitioning
		heom::partition_mapping partitioning;
		if (partition_file.empty()) {
			std::cout << "main(..): calculating partitioning..." << std::endl;
			// partitioning = create_simple_partitioning(complete_graph, rank, ranks);
			// partitioning = create_round_robin_partitioning(complete_graph, rank, ranks);
			partitioning = create_layer_partitioning(complete_graph, 0, ranks);
		} else {
			std::cout << "main(..): creating partitioning from file '" << partition_file << "'..." << std::endl;
			partitioning = create_file_partitioning(complete_graph, 0, ranks, partition_file);
		}

		for (rank_t r = 0; r < ranks; ++r) {
			std::cout << "main(..): creating hierarchy_partition for rank " << r << "..." << std::endl;
			heom::hierarchy_partition(complete_graph, partitioning, r, ranks);
			// TODO: add some nicer ouput here, and the needed members to heom::hierarchy_partition
			// TODO: collect information
		}
		// TODO: generate summary output with min/max, overall communication volume, cut edges, etc.

	} catch (...) {
		eptr = std::current_exception();
	}

	return heom::handle_main_exception(eptr);
}
