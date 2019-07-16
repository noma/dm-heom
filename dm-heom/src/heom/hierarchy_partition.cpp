// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/hierarchy_partition.hpp"

#include <set>

#include "debug.hpp"
#include "heom/config_error.hpp"

namespace heom {

// first partition has all layers and some of the last, i.e. first partition has all nodes as neighbours
partition_mapping create_simple_partitioning(const hierarchy_graph& graph, const rank_t& rank, const rank_t ranks)
{
	partition_mapping mapping(graph.nodes());

	rank_t nodes_per_rank = graph.nodes() / ranks;
	rank_t remainder = graph.nodes() % ranks;

	rank_t current_pos = 0;
	for (rank_t r = 0; r < ranks; ++r)
	{
		rank_t nodes_per_r = nodes_per_rank;
		// evenly distribute the remainder to the first remainder many ranks
		if (remainder > 0)
		{
			++nodes_per_r;
			--remainder;
		}

		// set rank for the next nodes_per_r nodes
		for (rank_t i = 0; i < nodes_per_r; ++i)
			mapping[current_pos + i] = r; // assign rank to each node

		current_pos += nodes_per_r;
	}

	return mapping;
}

// cycle through ranks
partition_mapping create_round_robin_partitioning(const hierarchy_graph& graph, const rank_t& rank, const rank_t ranks)
{
	partition_mapping mapping(graph.nodes());

	for (uid_t uid = 0; uid < graph.nodes(); ++uid)
		mapping[uid] = uid % ranks;

	return mapping;
}

// divide layers by ranks
partition_mapping create_layer_partitioning(const hierarchy_graph& graph, const rank_t& rank, const rank_t ranks)
{
	partition_mapping mapping(graph.nodes());

	mapping[0] = 0; // mandatory, 0th layer

	rank_t current_pos = 1;
	for (int_t current_depth = 1; current_depth <= graph.depth(); ++current_depth) {
		int_t nodes_per_current_depth = graph.nodes_per_depth()[current_depth];
		//for (rank_t i = 0; i < graph.nodes_per_depth()[current_depth]; ++i) {
		rank_t nodes_per_rank = nodes_per_current_depth / ranks;
		rank_t remainder = nodes_per_current_depth % ranks;

		for (rank_t r = 0; r < ranks; ++r)
		{
			rank_t nodes_per_r = nodes_per_rank;
			// evenly distribute the remainder to the first remainder many ranks
			if (remainder > 0)
			{
				++nodes_per_r;
				--remainder;
			}

			// set rank for the next nodes_per_r nodes
			for (rank_t i = 0; i < nodes_per_r; ++i)
				mapping[current_pos + i] = r; // assign rank to each node

			current_pos += nodes_per_r;
		}
	}

	return mapping;
}

partition_mapping create_file_partitioning(const hierarchy_graph& graph, const rank_t& rank, const rank_t ranks, const std::string& filename) {
	partition_mapping mapping;
	mapping.reserve(graph.nodes()); // prevent re-allocations

	std::ifstream file(filename);
	if (file.fail())
		throw config_error("could not open partition file: " + filename);

	rank_t max_rank = 0;
	for (std::string line; getline(file, line);) {
		rank_t r = boost::lexical_cast<rank_t>(line);
		mapping.push_back(r);
		max_rank = std::max(r, max_rank);
	}

	// check some basic properties of the loaded partitioning
	std::stringstream msg;
	if (mapping.size() != static_cast<size_t>(graph.nodes()))
		msg << "Partition mapping size (" << mapping.size() << ") does not match graph node count (" << graph.nodes() << ").";
	else if (max_rank > (ranks -1))
		msg << "Partition mapping's max rank (" << max_rank << ") exceeds rank count (" << ranks << ").";

	if (msg.str().size() > 0)
		throw std::runtime_error(msg.str());

	return mapping;
}



hierarchy_partition::hierarchy_partition(const hierarchy_graph& graph, const partition_mapping& mapping, const rank_t& rank, const rank_t ranks)
	: hierarchy_graph(graph.width(), graph.depth(), 0, 0, invalid_uid),
	  rank_(rank), ranks_(ranks), mapping_(mapping), rank_to_in_neighbours_(ranks), rank_to_out_neighbours_(ranks)
{
	// check edges between partitions and build data structures for the partitioning graph, and all neighbour relations
	int64_t in_neighbour_nodes = 0;
	int64_t out_neighbour_nodes = 0; // NOTE: included in compute_nodes_

	std::vector<uid_t> all_edges; // storage for combined plus and minus edges, re-used to avoid re-allocations
	all_edges.reserve(2 * graph.width());
	for (uid_t uid = 0; uid <  graph.nodes(); ++uid)
	{
		const rank_t uid_rank = mapping_[uid]; // rank of current node

		// create list of all edges
		all_edges.clear();
		const auto& plus_edges = graph.plus_edges()[uid];
		const auto& minus_edges = graph.minus_edges()[uid];
		all_edges.insert(all_edges.end(), plus_edges.begin(), plus_edges.end());
		all_edges.insert(all_edges.end(), minus_edges.begin(), minus_edges.end());

		// count number of computed nodes for this partition
		if (uid_rank == rank) // current node belongs to this partition
			++compute_nodes_;

		// iterate through adjacent nodes
		for (const auto& adj_uid : all_edges)
		{
			if (adj_uid != invalid_uid) // actual edge
			{
				rank_t adj_uid_rank = mapping_[adj_uid]; // partition rank for adj_uid

				if (uid_rank == rank) // current node belongs to this partition
				{
					if (adj_uid_rank != rank)  // adjacent node is not inside this partition
					{
						// collect for in-neighbours
						bool inserted = rank_to_in_neighbours_[adj_uid_rank].insert(adj_uid).second; // add remote adj_uid to neighbours needed *from* remote rank adj_uid for computation
						if (inserted)
							++in_neighbour_nodes; // keep count
						in_neighbour_ranks_.insert(adj_uid_rank); // build set of in-neighbour ranks for partition graph
					}

				}
				else // current node belongs to another partition
				{
					if (adj_uid_rank == rank) // adjacent uid is inside this partition
					{
						// collect for out-neighbours
						bool inserted = rank_to_out_neighbours_[uid_rank].insert(adj_uid).second; // add local adj_uid to neighbours needed by remote uid_rank
						if (inserted)
							++out_neighbour_nodes; // keep count
						out_neighbour_ranks_.insert(uid_rank);  // build set of out-neighbour ranks for partition graph
					}
				}
			}
		} // for all_edges
	} // for uids

	DEBUG_ONLY(
		std::cout << "rank " << rank << ": compute_nodes = " << compute_nodes_ << std::endl;
		std::cout << "rank " << rank << ": halo_nodes = " << in_neighbour_nodes << std::endl;
		std::cout << "rank " << rank << ": in_neighbour_ranks_.size() = " << in_neighbour_ranks_.size() << std::endl;
		std::cout << "rank " << rank << ": out_neighbour_ranks_.size() = " << out_neighbour_ranks_.size() << std::endl;
		for (const auto& pair : rank_to_in_neighbours_)
			std::cout << "rank " << rank << ": rank_to_in_neighbours_ for rank " << pair.first << " has size " << pair.second.size() << std::endl;
		for (const auto& pair : rank_to_out_neighbours_)
			std::cout << "rank " << rank << ": rank_to_out_neighbours_ for rank " << pair.first << " has size " << pair.second.size() << std::endl;
	)

	// TODO: the following section could be refactored into graph.subgraph(nodelist) and re-used for filtering
	// build up sub-graph

	// copy all nodes assigned to this rank, they are ordered by uid

	uid_to_lid_.insert(std::make_pair(invalid_uid, invalid_uid)); // invalid maps to invalid
	tuples_.reserve(compute_nodes_ + in_neighbour_nodes); // allocate memory to avoid re-allocations
	lid_to_uid_.reserve(compute_nodes_ + in_neighbour_nodes); // allocate memory to avoid re-allocations
	uid_t lid_count = 0; // count up local node ids for this partition

	for (uid_t uid = 0; uid <  graph.nodes(); ++uid)
	{
		if (mapping_[uid] == rank) // uid belongs to this partition
		{
			tuples_.push_back(graph.tuples()[uid]);
			lid_to_uid_.push_back(uid); // build up lok-up-table to get back to global id's from graph
			ASSERT_ONLY( bool inserted = )
				uid_to_lid_.insert(std::make_pair(uid, lid_count))
			ASSERT_ONLY( .second )
				;
			assert(inserted); // must be unique

			// set first_non_plus_uid
			// NOTE: assumption: local ids are still in global uid order
			if (uid >= graph.first_non_plus_uid() && first_non_plus_uid_ == invalid_uid) // encounter first uid without plus edges
				first_non_plus_uid_ = lid_count; // set this partition's first_non_plus_uid, which is a local id of the graph partition

			++lid_count;
			assert(lid_to_uid_.size() == static_cast<decltype(lid_to_uid_)::size_type>(lid_count));
			++nodes_;
		}
	}
	if (first_non_plus_uid_ == invalid_uid) // no uid with no plus edges in partition found
		first_non_plus_uid_ = lid_count; // == last assigned lid + 1

	// copy in-neighbour nodes from each neighbour rank, ordered by rank, ordered by uid per rank
	for (rank_t r = 0; r < ranks; ++r)
	{
		if(rank_to_in_neighbours_.find(r) != rank_to_in_neighbours_.end()) // check if rank is in-neighbour
		{
			ASSERT_ONLY( uid_t last_uid = invalid_uid; )
			for (const auto& adj_uid : rank_to_in_neighbours_.at(r)) // iterate over all neighbours from rank r
			{
				assert(adj_uid > last_uid); // make sure order is as expected
				tuples_.push_back(graph.tuples()[adj_uid]);
				lid_to_uid_.push_back(adj_uid);
				ASSERT_ONLY( bool inserted = )
					uid_to_lid_.insert(std::make_pair(adj_uid, lid_count))
				ASSERT_ONLY( .second )
					;
				assert(inserted); // must be unique

				++lid_count;
				assert(lid_to_uid_.size() == static_cast<decltype(lid_to_uid_)::size_type>(lid_count));
				++nodes_;

				ASSERT_ONLY( last_uid = adj_uid; )
			}
		}
	}

	assert(static_cast<decltype(nodes_)>(lid_count) == nodes_); // must be the same, we use two counters for readability and different types

	DEBUG_ONLY( std::cerr << "hierarchy_partition::hierarchy_partition(..): rank " << rank_ << "/" << ranks_
		                  << ": nodes = " << nodes_
	                      << ", compute_nodes_ = " << compute_nodes_
	                      << ", halo_nodes = " << in_neighbour_nodes << std::endl; )
	assert(nodes_ == compute_nodes_ + in_neighbour_nodes);
	assert(nodes() == halo_nodes() + compute_nodes());

	// copy edges and translate uid into local id (done separately as all lids for this partition's nodes are needed)

	plus_edges_.resize(nodes_, adj_list(width_));
	minus_edges_.resize(nodes_, adj_list(width_));

	for (uid_t lid = 0; lid < compute_nodes_; ++lid)
	{
		// corresponding adjacency list for current local id (lid) in original graph (uid)
		const auto& plus_adj_list = graph.plus_edges()[lid_to_uid_.at(lid)];
		const auto& minus_adj_list = graph.minus_edges()[lid_to_uid_.at(lid)];

		for (int64_t i = 0; i < width_; ++i)
		{
			// NOTE: use of at() to provoke exception instead of implicit insertion if not found
			plus_edges_[lid][i] = uid_to_lid_.at(plus_adj_list[i]);
			minus_edges_[lid][i] = uid_to_lid_.at(minus_adj_list[i]);
		}
	}

	// invalidate all edges for the neighbour nodes, as they are not used and most likely point outside the partition
	for (uid_t lid = compute_nodes_; lid < nodes_; ++lid)
	{
		for (int64_t i = 0; i < width_; ++i)
		{
			plus_edges_[lid][i] = invalid_uid;
			minus_edges_[lid][i] = invalid_uid;
		}

	}
}


} // namespace heom
