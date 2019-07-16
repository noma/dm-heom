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

int main(int argc, char* argv[])
{
	std::string heom_config_file = "population_dynamics.cfg";
	if (argc >= 2)
		heom_config_file = argv[1];

	std::string graph_file = heom_config_file + ".graph";
	if (argc >= 3)
		graph_file = argv[2];

	std::exception_ptr eptr;
	try {
		heom::population_dynamics_config heom_config(heom_config_file);

		// create hierarchy graph and partitioning
		heom::hierarchy_graph complete_graph(heom_config.baths_number(), heom_config.baths_matsubaras(), heom_config.system_ado_depth());

		std::cout << "Hierarchy graph stats:" << '\n'
		          << "\twidth:       " << complete_graph.width() << '\n'
		          << "\tdepth:       " << complete_graph.depth() << '\n'
		          << "\tnodes:       " << complete_graph.nodes() << '\n'
		          << "\tedges:       " << complete_graph.edges() << '\n'
		          << "\tplus_edges:  " << complete_graph.plus_edge_count() << '\n'
		          << "\tminus_edges: " << complete_graph.minus_edge_count() << std::endl;

		complete_graph.write_graph_file(graph_file);
		std::cout << "Graph file has been written to '" << graph_file << "'." << std::endl;
	} catch (...) {
		eptr = std::current_exception();
	}

	return heom::handle_main_exception(eptr);
}
