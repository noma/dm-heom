// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/hierarchy_graph.hpp"

#include <numeric>

#include <noma/debug/debug.hpp>

#include "heom/integer_partition.hpp"


namespace heom {

std::size_t tuple_hash(const tuple& v)
{
	return boost::hash_range(v.begin(), v.end());
}

const uid_t hierarchy_graph::invalid_uid = -2; // -2 for Mathematica compatibility

hierarchy_graph::hierarchy_graph(std::int64_t baths_number, std::int64_t baths_matsubaras, std::int64_t depth)
	: hierarchy_graph(baths_number * baths_matsubaras, depth)
{}

hierarchy_graph::hierarchy_graph(std::int64_t width, std::int64_t depth)
	: width_(width), depth_(depth)
{
	// binomial coefficient of (width + depth over ado_width), rounded to nearest integer, this is a size
	nodes_ = static_cast<decltype(nodes_)>(std::llround(boost::math::binomial_coefficient<double>(width_ + depth_, width_)));
	DEBUG_ONLY( std::cout << "hierarchy_graph::hierarchy_graph(..): nodes_ = " << nodes_ << std::endl; )
	compute_nodes_ = nodes_;

	// same as above, but without the last layer, which has no plus uids (i.e. downward edges), this is an uid, i.e. range from 0 to (nodes_ - 1)
	//first_non_plus_uid_ = static_cast<uid_t>(std::llround(boost::math::binomial_coefficient<double>(width_ + depth_ - 1, depth_ - 1))); // NOTE:  mathematically equivalent, but harder to read version
	first_non_plus_uid_ = static_cast<uid_t>(std::llround(boost::math::binomial_coefficient<double>(width_ + depth_ - 1, width_))); // first uid of last layer
	DEBUG_ONLY( std::cout << "hierarchy_graph::hierarchy_graph(..): first_non_plus_uid_ = " << first_non_plus_uid_ << std::endl; )
	assert(first_non_plus_uid_ > 0);


	// Step 1: Generate all nodes, as list of tuples (excitation vectors)
	tuples_.reserve(nodes_); // allocate memory
	nodes_per_depth_.resize(depth_ + 1, 0);

	tuple_to_uid_map tuple_to_uid(nodes_, tuple_hash); // excitation vector to uid map, the uid is counted up starting with 0 and establishes a simple numbering/indexing scheme

	int_t* int_part = nullptr;
	int_t int_part_size = 0;

	uid_t uid_counter = 0;

	// iterate over hierarchy levels top to bottom
	for (std::int64_t current_depth = 0; current_depth <= depth_; ++current_depth)
	{
		DEBUG_ONLY( std::cout << "hierarchy_graph::hierarchy_graph(..): constructing hierarchy nodes, level " << current_depth << " of " << depth_ << std::endl; )

		// fill in all integer partitions
		if (current_depth == 0) // special case: the zero result: fill in width_ 0's
		{
			int_part_size = 1;
			int_part = new int_t[width_];
			for (int_t j = 0; j < width_; ++j)
				int_part[j] = 0;
		}
		else
		{
			// compute all integer partitions for the current depth (physics: all ways to split n excitations on the baths)
			int_part_size = ip::size_my_integer_partitions(current_depth, width_);
			int_part = new int_t[int_part_size * width_];
			ip::my_integer_partitions(current_depth, width_, int_part);
		}

		// iterate over all integer partitions
#ifdef HEOM_MATHEMATICA_COMPAT
		// for Mathematica compatibility
		for (int_t i = int_part_size - 1; i >= 0; --i)
#else
		for (int_t i = 0; i < int_part_size; ++i)
#endif
		{
			// construct a tuple from an integer partition
			tuple tuple_ip;
			tuple_ip.resize(width_);
			for (std::int64_t j = 0; j < width_; ++j)
				tuple_ip[j] = int_part[i * width_ + j]; // element-wise copy

			// sort integer partition tuple
			std::sort(tuple_ip.begin(), tuple_ip.end());

			// iterate over all permutations of the sorted integer partition tuple
			do
			{
#ifdef HEOM_MATHEMATICA_COMPAT
				// reverse the key for better Mathematica compability
				// NOTE: we need a copy, as tuple_ip is systematically changed by std::next_permutation
				tuple tuple_copy(tuple_ip);
				std::reverse(tuple_copy.begin(), tuple_copy.end());
#else
				const tuple& tuple_copy = tuple_ip; // just create a const reference
#endif

				// fill tuples_, i.e. a tuple-valued list of vertices
				// same as: tuples_[uid_counter] = tuple_copy
				tuples_.push_back(tuple_copy);
				nodes_per_depth_[current_depth] += 1;

				// fill look-up table
				tuple_to_uid.insert(std::pair<tuple, int_t>(tuple_copy, uid_counter));

				uid_counter++;
				assert(uid_counter <= nodes_); // make sure nothing went wrong
			} while (std::next_permutation(tuple_ip.begin(), tuple_ip.end()));
		}
		delete[] int_part;
	}
	assert(uid_counter == nodes_);

	// Step 2: generate edges, i.e. the plus and minus indices

	// initialise with correct size of inner and outer list
	plus_edges_.resize(nodes_, adj_list(width_));
	minus_edges_.resize(nodes_, adj_list(width_));

	tuple plus_tuple(width_);
	tuple minus_tuple(width_);

	int_t plus_uid;
	int_t minus_uid;

	DEBUG_ONLY( std::cout << "hierarchy_graph::hierarchy_graph(..): creating plus/minus edges" << std::endl; )
	// generate plus/minus connections/indices for every node
	for (std::int64_t uid = 0; uid < nodes_; ++uid)
	{
		// read tuple for current matrix/node and compute its depth
		const tuple& current_tuple = tuples_[uid];
		int_t current_depth = hierarchy_graph::depth(current_tuple);

		// iterate over unit vectors (j-th element is 1)
		// and add/substract them from the current tuple to generate the connected plus/minus tuples
		for (std::int64_t j = 0; j < width_; ++j)
		{
			plus_tuple = current_tuple;
			minus_tuple = current_tuple;
			plus_tuple[j] += 1;
			minus_tuple[j] -= 1;

			if (current_depth < depth_) { // PLUS stays one level below the maximum depth, i.e. no downward edges in the last layer
				plus_uid = tuple_to_uid.at(plus_tuple); // look up the id for the resulting tuple
				++plus_edge_count_;
			} else {
				plus_uid = invalid_uid;
			}

			if (minus_tuple[j] < 0) { // MINUS cannot produce negative, i.e. invalid, index
				minus_uid = invalid_uid;
			} else {
				minus_uid = tuple_to_uid.at(minus_tuple);
				++minus_edge_count_;
			}

			plus_edges_[uid][j] = plus_uid;
			minus_edges_[uid][j] = minus_uid;
		}
	}

	// check first_non_plus_id for sanity: sum of tuple of uid = (first_non_plus - 1) of the complete graph must be exactly one less than those of the uid's one layer higher up
	DEBUG_ONLY(
		// access tuples
		auto first_non_plus_tuple = tuples()[first_non_plus_uid()];
		auto prior_tuple = tuples()[first_non_plus_uid() - 1]; // last tuple of layer above
		// sum over tuples
		auto first_non_plus_tuple_sum = std::accumulate(first_non_plus_tuple.begin(), first_non_plus_tuple.end(), 0);
		auto prior_tuple_sum = std::accumulate(prior_tuple.begin(), prior_tuple.end(), 0);
		// compare sums
		assert((first_non_plus_tuple_sum - 1) == prior_tuple_sum);
	)
}

hierarchy_graph::hierarchy_graph(int64_t width, int64_t depth, int64_t nodes, int64_t compute_nodes, uid_t first_non_plus_uid)
	: width_(width), depth_(depth), nodes_(nodes), compute_nodes_(compute_nodes), first_non_plus_uid_(first_non_plus_uid)
{
}

int_t hierarchy_graph::depth(const tuple& t) const
{
	return std::accumulate(t.begin(), t.end(), 0);
}

int_t hierarchy_graph::depth(const uid_t& uid) const
{
	const tuple& t = tuples_[uid];
	return depth(t);
}

void hierarchy_graph::write_graph_file(const std::string& filename)
{
	std::ofstream file(filename);

	if (!file.is_open())
		throw std::runtime_error("Could not open graph file for writing: " + filename);

	// generate body part, and count edges (as undirected, and only once for (u,v) and (v,u)) for the header
	size_t edges = 0;
	std::stringstream body;

	// build undirected graph data structure (vector of sets of adjacent nodes)
	std::vector<std::set<uid_t>> undirected_graph;
	undirected_graph.resize(nodes());

	for (uid_t uid = 0; uid < nodes(); ++uid) {
		// build adjacency list for undirected graph
		auto add_edges = [&] (const adj_list& adj) {
			for (const auto& adj_uid : adj) {
				if (adj_uid != invalid_uid) {
					undirected_graph[uid].insert(adj_uid);
					undirected_graph[adj_uid].insert(uid);
				}
			}
		}; // lambda add_edges

		add_edges(plus_edges()[uid]);
		add_edges(minus_edges()[uid]);
	}

	// count edges
	for (uid_t uid = 0; uid < nodes(); ++uid)
		edges += undirected_graph[uid].size();
	assert(edges % 2 == 0);
	edges /= 2;

	// write header
	file << nodes() << ' ' << edges;
	for (const auto& adj_set : undirected_graph) {
		file << '\n';
		for(const auto& adj_uid : adj_set) {
			file << ' ' << (adj_uid + 1); // NOTE: vertex numbering starts with 1 in graph file
		}
	}

}

} // namespace heom
