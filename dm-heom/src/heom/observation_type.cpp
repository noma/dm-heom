// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/observation_type.hpp"

#include <string>

#include <noma/typa/parser_error.hpp>

#include "heom/common.hpp"

namespace heom {

std::ostream& operator<<(std::ostream& out, const observation_type_t& t)
{
	out << observation_type_names.at(t);
	return out;
}

std::istream& operator>>(std::istream& in, observation_type_t& t)
{
	std::string value;
	std::getline(in, value);

	// get key to value	
	// NOTE: we trust observation_type_names to be complete here
	bool found = false;
	for (auto it = observation_type_names.begin(); it != observation_type_names.end(); ++it)
		if (it->second == value) {
			t = it->first;
			found = true;
			break;
		}

	if (!found)
		throw parser::parser_error("'" + value + "' is not a valid observation_type, allowed values: " + parser::get_map_value_list_string(observation_type_names));

	return in;
}

bool is_single_propagation_observation(const observation_type_t& obs_type)
{
	return obs_type == observation_type_t::matrix_diagonal ||
	       obs_type == observation_type_t::matrix_complete ||
	       obs_type == observation_type_t::matrix_trace_id ||
	       obs_type == observation_type_t::hierarchy_norm  ||
	       obs_type == observation_type_t::hierarchy_complete;
}


} // namespace heom

namespace noma {
namespace typa {

// see header for explanation
const std::string& observation_type_literal()
{
	auto gen_name_exp_str = []() {
		std::string name_exp_str;
		for (auto& name : heom::observation_type_names) {
			if (!name_exp_str.empty())
				name_exp_str += "|";
			name_exp_str += "\\b" + name.second + "\\b"; // NOTE: "\\b" is an '\b' regex assertion for a word boundary
		}
		name_exp_str = "(?:" + name_exp_str + ")";
		return name_exp_str;
	};

	static const std::string value { gen_name_exp_str() };

	return value;
}

} // namespace typa
} // namespace noma
