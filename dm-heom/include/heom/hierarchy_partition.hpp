// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_hierarchy_partition_hpp
#define heom_hierarchy_partition_hpp

#include "heom/hierarchy_graph.hpp"

#include <vector>

namespace heom {

using partition_mapping = std::vector<rank_t>;

partition_mapping create_simple_partitioning(const hierarchy_graph& graph, const rank_t& rank, const rank_t ranks);
partition_mapping create_round_robin_partitioning(const hierarchy_graph& graph, const rank_t& rank, const rank_t ranks);
partition_mapping create_layer_partitioning(const hierarchy_graph& graph, const rank_t& rank, const rank_t ranks);

/**
 * Read partition mapping from file and make sure it fits the graph and rank setup of the current run.
 *
 * @param graph the graph to partition.
 * @param rank rank of the calling entitie (e.g. MPI process).
 * @param ranks total number of ranks, i.e. partitions.
 * @param filename name/path of input file.
 * @return the partition mapping, i.e. a list of ranks, where the i-th value denotes the rank of the i-th node of the graph.
 */
partition_mapping create_file_partitioning(const hierarchy_graph& graph, const rank_t& rank, const rank_t ranks, const std::string& filename);

class hierarchy_partition : public hierarchy_graph
{
public:
	using adj_set = std::set<uid_t>;
	using neighbour_map = std::unordered_map<rank_t, adj_set>;

	hierarchy_partition(const hierarchy_graph& graph, const partition_mapping& mapping, const rank_t& rank, const rank_t ranks);

	rank_t rank() const { return rank_; }
	rank_t ranks() const { return ranks_; }

	const partition_mapping& mapping() const { return mapping_; }

	const std::vector<uid_t>& lid_to_uid() const { return lid_to_uid_; }
	const std::map<uid_t, uid_t>& uid_to_lid() const { return uid_to_lid_; }

	const std::set<rank_t>& in_neighbour_ranks() const { return in_neighbour_ranks_; }
	const std::set<rank_t>& out_neighbour_ranks() const { return out_neighbour_ranks_; }

	const neighbour_map& rank_to_in_neighbours() const { return rank_to_in_neighbours_; }
	const neighbour_map& rank_to_out_neighbours() const { return rank_to_out_neighbours_; }

private:
	rank_t rank_;
	rank_t ranks_;

	partition_mapping mapping_;

	std::vector<uid_t> lid_to_uid_; // keeps the link from the partitions sequential node id's to the graphs vertex id's it was constructed from
	std::map<uid_t, uid_t> uid_to_lid_;

	std::set<rank_t> in_neighbour_ranks_; // ranks that have vertices this partition needs for computation
	std::set<rank_t> out_neighbour_ranks_; // ranks that need vertices from this partition for computation

	// get list of hierarchy_nodes to each neighbour rank, for in/out
	neighbour_map rank_to_in_neighbours_;
	neighbour_map rank_to_out_neighbours_;
};

} // namesapce heom

#endif // heom_hierarchy_partition_hpp
