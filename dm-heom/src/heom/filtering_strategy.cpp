// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
// Copyright (c) 2017 Lisa Gaedke-Merzhäuser, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/filtering_strategy.hpp"

#include <string>

#include <noma/typa/parser_error.hpp>
#include <noma/typa/util.hpp>

#include "heom/common.hpp"

namespace heom {

std::ostream& operator<<(std::ostream& out, const filtering_strategy_t& s)
{
	out << filtering_strategy_names.at(s);
	return out;
}

std::istream& operator>>(std::istream& in, filtering_strategy_t& s)
{
	std::string value;
	std::getline(in, value);

	// get key to value
	// NOTE: we trust filtering_strategy_names to be complete here
	bool found = false;
	// NOTE: iterate through list of pairs <filtering_strategy_t, std::string>, denoted by <first,second>
	for (auto it = filtering_strategy_names.begin(); it != filtering_strategy_names.end(); ++it)
		if (it->second == value) {
			s = it->first;
			found = true;
			break;
		}

	if (!found)
		throw parser::parser_error("'" + value + "' is not a valid filtering strategy, allowed strategies: " + parser::get_map_value_list_string(filtering_strategy_names));

	return in;
}

} // namespace heom