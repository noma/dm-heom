// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
// Copyright (c) 2017 Lisa Gaedke-Merzhäuser, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_hierarchy_mask_hpp
#define heom_hierarchy_mask_hpp

#include <vector>
#include "heom/hierarchy_graph.hpp"
#include "heom/filtering_strategy.hpp"

namespace heom {

enum class mask_t : uint8_t {
	MASK_FALSE = 0,
	MASK_TRUE = 1,
	MASK_PROTECTED = 2,
};

using hierarchy_mask = std::vector<mask_t>;

// define switch function to go to the correct create filter function
hierarchy_mask create_hierarchy_mask(filtering_strategy_t filtering_strategy, const hierarchy_graph& complete_graph, int_t bath_masubaras, int_t filtering_first_layer);

// if filtering strategy == NONE
// NOTE: second argument is unused
hierarchy_mask create_complete_hierarchy_mask(const hierarchy_graph& complete_graph, int_t, int_t);

// otherwise
hierarchy_mask create_single_excitation_mask(const hierarchy_graph& complete_graph, int_t bath_masubaras, int_t filtering_first_layer);

// outputs hierarchy mask buffer to stream
void write_hierarchy_mask(const mask_t* host_buffer, size_t size_array, std::ostream& stream);

// outputs hierarchy mask counter
void write_hierarchy_mask_stats(const hierarchy_mask& mask, std::ostream& os);

} // namespace heom

#endif // heom_hierarchy_mask_hpp
