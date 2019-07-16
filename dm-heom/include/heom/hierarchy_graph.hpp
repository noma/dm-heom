// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_hierarchy_graph_hpp
#define heom_hierarchy_graph_hpp

#include <unordered_map>

#include <boost/functional/hash.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "heom/common.hpp"


namespace heom {

using uid_t = std::int32_t; // unique node id

// represents a tuple in the HEOM sense (kind of an excitation vector) which uniquely describes one ado_tuple/matrix/hierarchy node
using tuple = std::vector<int_t>; // excitation vector, value type of each node
// hash function for tuple
//auto tuple_hash = [](const tuple& v) { return boost::hash_range(v.begin(), v.end()); }; // TODO: report compiler bug, causes multiple definition linker error with Intel compiler, tested with 16.0.3 and 17.0.1
std::size_t tuple_hash(const tuple& v);

class hierarchy_graph
{
public:
	//using tuple_to_uid_map = std::unordered_map<tuple, uid_t, decltype(tuple_hash)>;
	using tuple_to_uid_map = std::unordered_map<tuple, uid_t, std::function<decltype(tuple_hash)>>;
	using node_list = std::vector<tuple>;
	using adj_list = std::vector<uid_t>;
	using edge_list = std::vector<adj_list>;

	static const uid_t invalid_uid;

	hierarchy_graph(std::int64_t baths_number, std::int64_t baths_matsubaras, std::int64_t depth);
	hierarchy_graph(std::int64_t width, std::int64_t depth);

	int64_t width() const { return width_; }
	int64_t depth() const { return depth_; }
	int64_t nodes() const { return nodes_; }
	int64_t compute_nodes() const { return compute_nodes_; }
	int64_t halo_nodes() const { return nodes_ - compute_nodes_; }
	int64_t edges() const { return plus_edge_count_ + minus_edge_count_; }
	int64_t plus_edge_count() const { return plus_edge_count_; }
	int64_t minus_edge_count() const { return minus_edge_count_; }

	uid_t first_non_plus_uid() const { return first_non_plus_uid_; }

	const node_list& tuples() const { return tuples_; }
	const edge_list& plus_edges() const { return plus_edges_; }
	const edge_list& minus_edges() const { return minus_edges_; }

	const std::vector<int_t>& nodes_per_depth() const { return nodes_per_depth_; }

	// compute depth of a tuple
	int_t depth(const tuple& t) const;
	// compute depth of an uid
	int_t depth(const uid_t& uid) const;

	/**
	 * Write graph to file, METIS and Boost.Graph compatible adjacency list
	 * format, i.e. first row contains two space separated numbers for the
	 * node and edge count. Every other row i cantains the adjacency list for
	 * node i (count starting at 1).
	 *
	 * @param filename name/path of the output file.
	 */
	void write_graph_file(const std::string& filename);

protected:
	hierarchy_graph() = default;
	hierarchy_graph(int64_t width, int64_t depth, int64_t nodes, int64_t compute_nodes, uid_t first_non_plus_uid);

//private:
	int64_t width_ = 0; // max. number of outgoing plus or minus edges per node
	int64_t depth_ = 0; // number of hierarchy layers within the graph
	int64_t nodes_ = 0; // number of nodes
	int64_t compute_nodes_ = 0; // nodes to be computed
	int64_t plus_edge_count_ = 0; // number of plus edges
	int64_t minus_edge_count_ = 0; // number of minu edges

	// first node that has no outgoing downward (plus) edges
	uid_t first_non_plus_uid_ = invalid_uid;

	node_list tuples_;
	edge_list plus_edges_;
	edge_list minus_edges_;

	std::vector<int_t> nodes_per_depth_;
};

} // namespace heom

#endif // heom_hierarchy_graph_hpp
