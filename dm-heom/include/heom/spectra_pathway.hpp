// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef spectra_pathway_hpp
#define spectra_pathway_hpp

#include "heom/pathway.hpp"

namespace heom {

enum class spectra_pathway_t
{
	gbrp = 0,
	serp,
	esarp,
	gbnr,
	senr,
	esanr
};

// NOTE: names need to be unique and not a prefix of one another
const std::map<spectra_pathway_t, std::string> spectra_pathway_names {
	{ spectra_pathway_t::gbrp,  "gbrp"  },
	{ spectra_pathway_t::serp,  "serp"  },
	{ spectra_pathway_t::esarp, "esarp" },
	{ spectra_pathway_t::gbnr,  "gbnr"  },
	{ spectra_pathway_t::senr,  "senr"  },
	{ spectra_pathway_t::esanr, "esanr" }
};

using spectra_pathway_specification = pathway_specification<4>;

using spectra_pathway_specification_map = std::map<spectra_pathway_t, spectra_pathway_specification>;

/**
 * first set: +1/-1 = dipole_plus/dipole_minus, second set multiply from (R)ight or (L)eft
 * REPHASING paths,
 * GBRP  = +{{-, +, +, -}, {R, R, L, L}}
 * SERP  = +{{-, +, +, -}, {R, L, R, L}}
 * ESARP = -{{-, +, +, -}, {R, L, L, L}}
 * NONREPHASING paths
 * GBNR  = +{{+, -, +, -}, {L, L, L, L}}
 * SENR  = +{{+, -, +, -}, {L, R, R, L}}
 * ESANR = -{{+, -, +, -}, {L, R, L, L}}
 *
 * Every Pathway has a fourth entry for the trace, that is always -1 and Left
 */
const spectra_pathway_specification_map spectra_pathway_to_specification {
	{spectra_pathway_t::gbrp,  {{{{pathway_dipole_sign_t::minus, pathway_dipole_side_t::right},
	                              {pathway_dipole_sign_t::plus, pathway_dipole_side_t::right},
	                              {pathway_dipole_sign_t::plus, pathway_dipole_side_t::left},
	                              {pathway_dipole_sign_t::minus, pathway_dipole_side_t::left}}},
	                            pathway_sign_t::plus, sites_to_states_mode_t::with_ground_state}},

	{spectra_pathway_t::serp,  {{{{pathway_dipole_sign_t::minus, pathway_dipole_side_t::right},
	                              {pathway_dipole_sign_t::plus, pathway_dipole_side_t::left},
	                              {pathway_dipole_sign_t::plus, pathway_dipole_side_t::right},
	                              {pathway_dipole_sign_t::minus, pathway_dipole_side_t::left}}},
	                            pathway_sign_t::plus, sites_to_states_mode_t::with_ground_state}},

	{spectra_pathway_t::esarp, {{{{pathway_dipole_sign_t::minus, pathway_dipole_side_t::right},
	                              {pathway_dipole_sign_t::plus, pathway_dipole_side_t::left},
	                              {pathway_dipole_sign_t::plus, pathway_dipole_side_t::left},
	                              {pathway_dipole_sign_t::minus, pathway_dipole_side_t::left}}},
	                            pathway_sign_t::minus, sites_to_states_mode_t::excited_state_absorption}},

	{spectra_pathway_t::gbnr,  {{{{pathway_dipole_sign_t::plus, pathway_dipole_side_t::left},
	                              {pathway_dipole_sign_t::minus, pathway_dipole_side_t::left},
	                              {pathway_dipole_sign_t::plus, pathway_dipole_side_t::left},
	                              {pathway_dipole_sign_t::minus, pathway_dipole_side_t::left}}},
	                            pathway_sign_t::plus, sites_to_states_mode_t::with_ground_state}},

	{spectra_pathway_t::senr,  {{{{pathway_dipole_sign_t::plus, pathway_dipole_side_t::left},
	                              {pathway_dipole_sign_t::minus, pathway_dipole_side_t::right},
	                              {pathway_dipole_sign_t::plus, pathway_dipole_side_t::right},
	                              {pathway_dipole_sign_t::minus, pathway_dipole_side_t::left}}},
	                            pathway_sign_t::plus, sites_to_states_mode_t::with_ground_state}},

	{spectra_pathway_t::esanr, {{{{pathway_dipole_sign_t::plus, pathway_dipole_side_t::left},
	                              {pathway_dipole_sign_t::minus, pathway_dipole_side_t::right},
	                              {pathway_dipole_sign_t::plus, pathway_dipole_side_t::left},
	                              {pathway_dipole_sign_t::minus, pathway_dipole_side_t::left}}},
	                            pathway_sign_t::minus, sites_to_states_mode_t::excited_state_absorption}},
};

std::ostream& operator<<(std::ostream& out, const spectra_pathway_t& t);
std::istream& operator>>(std::istream& in, spectra_pathway_t& t);

} // namespace heom


namespace noma {
namespace typa {

const std::string& spectra_pathway_literal();

template<>
struct type_to_regexp<heom::spectra_pathway_t>
{
	static const std::string& exp_str() { return spectra_pathway_literal(); }
};


} // namespace typa
} // namespace noma


#endif // spectra_pathway_hpp
