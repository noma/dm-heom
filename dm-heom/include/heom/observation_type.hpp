// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_observation_type_hpp
#define heom_observation_type_hpp

#include <iostream>
#include <map>

#include <noma/typa/typa.hpp>

namespace heom {

enum class observation_type_t
{
	matrix_diagonal = 0,
	matrix_complete,
	matrix_complete_rate_kernel,
	matrix_trace_id,
	matrix_trace_linear_absorption,
	matrix_trace_transient_absorption,
	matrix_trace_static_fluorescence,
	matrix_trace_two_dimensional_spectra,
	hierarchy_norm,
	hierarchy_complete
};

// NOTE: names need to be unique and not a prefix of one another
const std::map<observation_type_t, std::string> observation_type_names {
	{ observation_type_t::matrix_diagonal, "matrix_diagonal" },
	{ observation_type_t::matrix_complete, "matrix_complete" },
	{ observation_type_t::matrix_complete_rate_kernel, "matrix_complete_rate_kernel" },
	{ observation_type_t::matrix_trace_id, "matrix_trace_id" },
	{ observation_type_t::matrix_trace_linear_absorption, "matrix_trace_linear_absorption" },
	{ observation_type_t::matrix_trace_transient_absorption, "matrix_trace_transient_absorption" },
	{ observation_type_t::matrix_trace_static_fluorescence, "matrix_trace_static_fluorescence" },
	{ observation_type_t::matrix_trace_two_dimensional_spectra, "matrix_trace_two_dimensional_spectra" },
	{ observation_type_t::hierarchy_norm, "hierarchy_norm" },
	{ observation_type_t::hierarchy_complete, "hierarchy_complete" }
};

std::ostream& operator<<(std::ostream& out, const observation_type_t& t);
std::istream& operator>>(std::istream& in, observation_type_t& t);

bool is_single_propagation_observation(const observation_type_t& obs_type);

} // namespace heom


namespace noma {
namespace typa {

const std::string& observation_type_literal();

template<>
struct type_to_regexp<heom::observation_type_t>
{
	static const std::string& exp_str() { return observation_type_literal(); }
};


} // namespace typa
} // namespace noma


#endif // heom_observation_type_hpp
