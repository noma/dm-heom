// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/program_task.hpp"

#include <string>

#include <noma/typa/parser_error.hpp>
#include <noma/typa/util.hpp>

#include "heom/common.hpp"

namespace heom {

std::ostream& operator<<(std::ostream& out, const program_task_t& t)
{
	out << program_task_names.at(t);
	return out;
}

std::istream& operator>>(std::istream& in, program_task_t& t)
{
	std::string value;
	std::getline(in, value);

	// get key to value	
	// NOTE: we trust program_task_names to be complete here
	bool found = false;
	for (auto it = program_task_names.begin(); it != program_task_names.end(); ++it)
		if (it->second == value) {
			t = it->first;
			found = true;
			break;
		}

	if (!found)
		throw parser::parser_error("'" + value + "' is not a valid program_task, allowed values: " + parser::get_map_value_list_string(program_task_names));

	return in;
}

} // namespace heom
