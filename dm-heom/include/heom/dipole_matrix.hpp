// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_dipole_matrix_hpp
#define heom_dipole_matrix_hpp

#include "heom/common.hpp"
#include "heom/dipole_config_entries.hpp"
#include "heom/laser_config_entries.hpp"
#include "heom/pathway.hpp"
#include "heom/sites_to_states.hpp"

namespace heom {

real_vector_t cross_product(const real_vector_t& a, const real_vector_t& b);

real_vector_t compute_dipole_vector(int_t sites,
                                    const heom::dipole_config_entries::dipole_directions_type& dipole_directions,
                                    const heom::laser_config_entries::laser_directions_type& laser_directions,
                                    int_t laser_direction_index,
                                    real_t polarization_angle);

complex_matrix_t compute_dipole_matrix(const real_vector_t& dipole_vector,
                                       sites_to_states_mode_t mode);

template<class PATHWAY_T, bool CIRCULAR = false>
complex_matrix_t get_dipole_matrix_for_pathway(const PATHWAY_T& pathway_spec, size_t spec_index, int_t system_sites, const dipole_config_entries& dipole_config, size_t tensor_index)
{
	auto sts_mode = pathway_spec.sites_to_states_mode();

	// current tensor component
	const int_t tensor_component = dipole_config.dipole_tensor_components().at(tensor_index).at(spec_index);

	heom::real_vector_t dipole_vector(system_sites);

	if (!CIRCULAR) {
		for (size_t i = 0; i < dipole_vector.size(); ++i)
			dipole_vector[i] = dipole_config.dipole_directions()[i][tensor_component] * dipole_config.dipole_strengths()[i];
	} else {
		// circular dichroism case m^+/- vs. µ^+/-
		for (size_t i = 0; i < dipole_vector.size(); ++i) {
			dipole_vector[i] = cross_product(dipole_config.dipole_centers()[i], dipole_config.dipole_directions()[i])[tensor_component] * dipole_config.dipole_strengths()[i];
		}
	}

	auto dipole_matrix = heom::compute_dipole_matrix(dipole_vector, sts_mode);

	// transpose if dipole_plus is required
	if (pathway_spec.pathway()[spec_index].first == heom::pathway_dipole_sign_t::plus)
		dipole_matrix = dipole_matrix.transposed();

	return dipole_matrix;
}


template<class PATHWAY_T>
complex_matrix_t get_dipole_matrix_for_pathway_circular(const PATHWAY_T& pathway_spec, size_t spec_index, int_t system_sites, const dipole_config_entries& dipole_config, size_t tensor_index)
{
	return get_dipole_matrix_for_pathway<PATHWAY_T, true>(pathway_spec, spec_index, system_sites, dipole_config, tensor_index);
}


template<class PATHWAY_T>
complex_t get_dipole_tensor_prefactor(const PATHWAY_T& pathway_spec, const dipole_config_entries& dipole_config, size_t tensor_index)
{
	return (pathway_spec.sign() == heom::pathway_sign_t::plus)
	      ?   dipole_config.dipole_tensor_prefactors().at(tensor_index)
	      : - dipole_config.dipole_tensor_prefactors().at(tensor_index);

}

} // namespace heom

#endif // heom_dipole_matrix_hpp