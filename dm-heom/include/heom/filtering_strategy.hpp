// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
// Copyright (c) 2017 Lisa Gaedke-Merzhäuser, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_filtering_strategy_hpp
#define heom_filtering_strategy_hpp

#include <iostream>
#include <map>

namespace heom {

enum class filtering_strategy_t
{
	NONE = 0,
	SINGLE_EXCITATION
};

const std::map<filtering_strategy_t, std::string> filtering_strategy_names {
	{ filtering_strategy_t::NONE, "none" },
	{ filtering_strategy_t::SINGLE_EXCITATION, "single_excitation" }};

std::ostream& operator<<(std::ostream& out, const filtering_strategy_t& s);
std::istream& operator>>(std::istream& in, filtering_strategy_t& s);

} // namespace heom

#endif // heom_filtering_strategy_hpp
