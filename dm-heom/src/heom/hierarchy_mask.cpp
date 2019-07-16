// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
// Copyright (c) 2017 Lisa Gaedke-Merzhäuser, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/hierarchy_mask.hpp"

#include <numeric>

namespace heom {

using hierarchy_mask = std::vector<mask_t>;

// define switch function to go to the correct create filter function

// int_t bath_masubaras,
hierarchy_mask create_hierarchy_mask(filtering_strategy_t filtering_strategy, const hierarchy_graph& complete_graph, int_t baths_matsubaras, int_t filtering_first_layer)
{
	switch (filtering_strategy) {
		case filtering_strategy_t::NONE:
			return create_complete_hierarchy_mask(complete_graph, baths_matsubaras, filtering_first_layer);
		case filtering_strategy_t::SINGLE_EXCITATION:
			return create_single_excitation_mask(complete_graph, baths_matsubaras, filtering_first_layer);
		default:
			throw std::runtime_error("create_hierarchy_mask(): error: unknown filtering strategy");
	}
}

// if filtering strategy == NONE
// int_t,
hierarchy_mask create_complete_hierarchy_mask(const hierarchy_graph& complete_graph, int_t, int_t)
{
	hierarchy_mask new_mask;
	new_mask.resize(complete_graph.nodes());
	for (auto& m : new_mask) {
		m = mask_t::MASK_TRUE;
	}
	return new_mask;
}

// otherwise
hierarchy_mask create_single_excitation_mask(const hierarchy_graph& complete_graph, int_t baths_matsubaras, int_t filtering_first_layer)
{
	// we have a single excitation whenever excitations belong to the same bath (matsubara property)
	hierarchy_mask new_mask;
	new_mask.resize(complete_graph.nodes());
	size_t mask_index = 0;
	size_t false_count = 0;
	size_t true_count = 0;

	// get length of a tuple,
	int_t tuple_length = complete_graph.width();
	int_t baths = tuple_length / baths_matsubaras;

	for (auto& tuple : complete_graph.tuples()) {
		// NOTE: tuple is std::vector<int>

		if (complete_graph.depth(tuple) < filtering_first_layer) {
			new_mask[mask_index] = mask_t::MASK_TRUE;
			++true_count;
			++mask_index;
		} else { // for higher layers initialise with MASK_FALSE and overwrite below, when we want MASK_TRUE.

			int_t elem_sum = std::accumulate(tuple.begin(), tuple.end(), 0);

			// consider tuple-entries of each bath
			bool found = false;
			for (int_t bath = 0; bath < baths; ++bath) {
				// sum over all entries in a tuple which belong to the same bath
				// NOTE: unfortunately this heavily relies on the particular ordering of the entries in the tuple
				int_t bath_sum = std::accumulate(tuple.begin() + bath * baths_matsubaras, tuple.begin() + (bath + 1) * baths_matsubaras, 0);

				if (bath_sum == elem_sum) {
					found = true;
					break;
				}
			}

			if (found) {
				new_mask[mask_index] = mask_t::MASK_TRUE;
				++true_count;
				DEBUG_ONLY(
					std::cout << "heom::create_single_excitation_mask(): new tuple: ";
					for (size_t i = 0; i < tuple.size(); ++i)
						std::cout << tuple[i] << " ";
				)
			} else {
				new_mask[mask_index] = mask_t::MASK_FALSE;
				++false_count;
			}

			++mask_index;
		}
	}

	DEBUG_ONLY(
	           std::cout << "heom::create_single_excitation_mask(): " << std::endl;
	           std::cout << default_delimiter << "tuple_length: " << tuple_length << std::endl;
	           std::cout << default_delimiter << "baths:        " << baths << std::endl;
	           std::cout << default_delimiter << "true count:   " << true_count << std::endl;
	           std::cout << default_delimiter << "false count:  " << false_count << std::endl;
	           std::cout << default_delimiter << "mask index:   " << mask_index << std::endl;
	)

	return new_mask;
}


void write_hierarchy_mask(const mask_t* host_buffer, size_t size_array, std::ostream& stream)
{
	auto format = int_format;
	for (size_t i = 0; i < size_array; ++i) {
		stream << format % static_cast<int>(host_buffer[i]);
		if (i != size_array - 1) { stream << default_delimiter; }
	}
	stream << std::endl;
}

// hierarchy mask stats
void write_hierarchy_mask_stats(const hierarchy_mask& mask, std::ostream& os)
{
	int_t true_count = 0;
	int_t protected_count = 0;
	int_t false_count =0;
	int_t total_number_of_matrices = mask.size();

	for (auto& m : mask) {
		switch (m) {
			case mask_t::MASK_TRUE:
				++true_count;
				break;
			case mask_t::MASK_PROTECTED:
				++protected_count;
				break;
			case mask_t::MASK_FALSE:
				++false_count;
				break;
			default:
				throw std::runtime_error("write_hierarchy_mask_stats(): error: encountered unknown hierarchy_mask type");
		}

	}

	real_t false_count_percent = (real_t) false_count / (real_t) total_number_of_matrices * 100;
	real_t true_count_percent = (real_t) true_count / (real_t) total_number_of_matrices * 100;
	real_t protected_count_percent = (real_t) protected_count / (real_t) total_number_of_matrices * 100;

	int_t total_active = true_count + protected_count;
	real_t total_active_percent = (real_t) total_active / (real_t) total_number_of_matrices * 100;

	size_t space_count_header = 17;
	std::string first_table_space_count = "%10i ";
	std::string second_table_space_count = "%22.2f";

	os << "matrices" << std::string(space_count_header, ' ') << "absolute" << std::string(space_count_header, ' ') << "relative" << std::endl;
	os << "------------------------------------------------------------" << std::endl;
	os << "total matrices:        "
	   << boost::format(first_table_space_count) % total_number_of_matrices
	   << boost::format(second_table_space_count) % 100.00 <<" %" << std::endl;
	os << "total active matrices: "
	   << boost::format(first_table_space_count) % total_active
	   << boost::format(second_table_space_count) % total_active_percent << " %" << std::endl;
	os << "mask false:            "
	   << boost::format(first_table_space_count) % false_count
	   << boost::format(second_table_space_count) % false_count_percent <<" %" << std::endl;
	os << "mask true:             "
	   << boost::format(first_table_space_count) % true_count
	   << boost::format(second_table_space_count) % true_count_percent <<" %" << std::endl;
	os << "mask protected:        "
	   << boost::format(first_table_space_count) % protected_count
	   << boost::format(second_table_space_count) % protected_count_percent <<" %" << std::endl;
}

} // namespace heom

