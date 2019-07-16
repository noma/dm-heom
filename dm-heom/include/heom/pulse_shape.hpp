// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_pulse_shape_hpp
#define heom_pulse_shape_hpp

#include <iostream>
#include <map>

namespace heom {

enum class pulse_shape_t {
	GAUSSIAN
};

const std::map<pulse_shape_t, std::string> pulse_shape_names {
	{ pulse_shape_t::GAUSSIAN, "Gaussian" }};

std::ostream& operator<<(std::ostream& out, const pulse_shape_t& ps);
std::istream& operator>>(std::istream& in, pulse_shape_t& ps);


} // namespace heom

#endif // heom_pulse_shape_hpp
