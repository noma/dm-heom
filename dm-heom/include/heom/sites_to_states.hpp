// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_sites_to_states_hpp
#define heom_sites_to_states_hpp

#include "heom/common.hpp"

namespace heom {

enum class sites_to_states_mode_t
{
	identity, // one exciton
	with_ground_state, // add ground state to identity, i.e. zero excitons
	excited_state_absorption // add ground state and states for two excitons
};


int_t sites_to_states(int_t sites, sites_to_states_mode_t mode);

// TODO: more sophisticated version that takes couplings into account and finds a minimal baths_per_state
int_t max_baths_per_state(int_t max_baths_per_site, sites_to_states_mode_t mode);

complex_matrix_t state_hamiltonian(const complex_matrix_t& hamiltonian, sites_to_states_mode_t mode);

parser::matrix<int_t> state_baths_coupling(const parser::matrix<int_t>& coupling, sites_to_states_mode_t mode);

} // namespace heom

#endif // heom_sites_to_states