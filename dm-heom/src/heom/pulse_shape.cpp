// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/pulse_shape.hpp"

#include <string>

#include <noma/typa/parser_error.hpp>
#include <noma/typa/util.hpp>

#include "heom/common.hpp"

namespace heom {

std::ostream& operator<<(std::ostream& out, const pulse_shape_t& t)
{
	out << pulse_shape_names.at(t);
	return out;
}

std::istream& operator>>(std::istream& in, pulse_shape_t& t)
{
	std::string value;
	std::getline(in, value);

	// get key to value
	// NOTE: we trust pulse_shape_names to be complete here
	bool found = false;
	for (auto it = pulse_shape_names.begin(); it != pulse_shape_names.end(); ++it)
		if (it->second == value) {
			t = it->first;
			found = true;
			break;
		}

	if (!found)
		throw parser::parser_error("'" + pulse_shape_names.at(t) + "' is not a valid pulse_shape.");

	return in;
}

} // namespace heom
