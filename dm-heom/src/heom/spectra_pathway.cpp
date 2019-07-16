// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/spectra_pathway.hpp"

#include <string>

#include <noma/typa/parser_error.hpp>

#include "heom/common.hpp"

namespace heom {

std::ostream& operator<<(std::ostream& out, const spectra_pathway_t& t)
{
	out << spectra_pathway_names.at(t);
	return out;
}

std::istream& operator>>(std::istream& in, spectra_pathway_t& t)
{
	std::string value;
	std::getline(in, value);

	// get key to value	
	// NOTE: we trust spectra_pathway_names to be complete here
	bool found = false;
	for (auto it = spectra_pathway_names.begin(); it != spectra_pathway_names.end(); ++it)
		if (it->second == value) {
			t = it->first;
			found = true;
			break;
		}

	if (!found)
		throw parser::parser_error("'" + value + "' is not a valid spectra_pathways, allowed values: " + parser::get_map_value_list_string(spectra_pathway_names));

	return in;
}

} // namespace heom

namespace noma {
namespace typa {

// see header for explanation
const std::string& spectra_pathway_literal()
{
	auto gen_name_exp_str = []() {
		std::string name_exp_str;
		for (auto& name : heom::spectra_pathway_names) {
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
