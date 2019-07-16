// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef dipole_pseudo_pathway_hpp
#define dipole_pseudo_pathway_hpp

#include "heom/dipole_config_entries.hpp"
#include "heom/pathway.hpp"

namespace heom {

/**
 * Generic pseudo pathway for dipole calculation = +{{-, +}, {L, L}}
 */

using dipole_pseudo_pathway_specification = pathway_specification<2>;

const dipole_pseudo_pathway_specification dipole_pseudo_pathway_wgs_spec {
	{{{pathway_dipole_sign_t::minus,  pathway_dipole_side_t::left},
	  {pathway_dipole_sign_t::plus, pathway_dipole_side_t::left}}},
	pathway_sign_t::plus,
	sites_to_states_mode_t::with_ground_state
};

const dipole_pseudo_pathway_specification dipole_pseudo_pathway_esa_spec {
	{{{pathway_dipole_sign_t::minus,  pathway_dipole_side_t::left},
	  {pathway_dipole_sign_t::plus, pathway_dipole_side_t::left}}},
	pathway_sign_t::plus,
	sites_to_states_mode_t::excited_state_absorption
};

} // namespace heom

#endif // dipole_pseudo_pathway_hpp
