// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/dipole_matrix.hpp"

#include "heom/constants.hpp"

namespace heom {

real_vector_t cross_product(const real_vector_t& a, const real_vector_t& b)
{
	assert(a.size() == 3 && b.size() == 3);

	real_vector_t c { 3 };

	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];

	return c;
}

real_t scalar_product(const real_vector_t& a, const real_vector_t& b)
{
	assert(a.size() == b.size());

	real_t c = 0.0;

	for (size_t i = 0; i < a.size(); ++i) {
		c += a[i] * b[i];
	}

	return c;
}

real_vector_t operator+(const real_vector_t& a, const real_vector_t& b)
{
	assert(a.size() == b.size());

	real_vector_t c { a.size(), 0.0 };

	for (size_t i = 0; i < a.size(); ++i) {
		c[i] = a[i] + b[i];
	}

	return c;
}

real_vector_t operator*(real_t s, const real_vector_t& a)
{
	real_vector_t c { a.size(), 0.0 };

	for (size_t i = 0; i < a.size(); ++i) {
		c[i] = s * a[i];
	}

	return c;
}

real_vector_t norm(const real_vector_t a)
{
	real_vector_t n;

	return (1.0 / sqrt(scalar_product(a, a))) * a;
}

real_vector_t compute_dipole_vector(int_t sites,
                                    const heom::dipole_config_entries::dipole_directions_type& dipole_directions,
                                    const heom::laser_config_entries::laser_directions_type& laser_directions,
                                    int_t laser_direction_index, // which laser direction to use
                                    real_t polarization_angle) // overall rotation angle
{
	real_vector_t dipole_vector(sites); // µS in Mathematica prototype

	polarization_angle = constants::pi * polarization_angle / 180.0; // convert to rad

	// compute e-field (cartesian direction vector)
	const real_vector_t& laser_direction = laser_directions[laser_direction_index];
	const real_vector_t& laser_direction_other = laser_directions[(laser_direction_index + 1) % laser_directions.size()];

	// base vector, that is perpendicular to the current and another laser direction
	const real_vector_t e_field_base = norm(cross_product(laser_direction, laser_direction_other));

	const real_vector_t laser_direction_normed = norm(laser_direction);
	// rotate e_field_base by polarization_angle around laser_direction
	const real_vector_t e_field = scalar_product(laser_direction_normed, e_field_base) * laser_direction_normed
	                            + cos(polarization_angle) * cross_product(cross_product(laser_direction_normed, e_field_base), laser_direction_normed)
	                            + sin(polarization_angle) * cross_product(laser_direction_normed, e_field_base);

	DEBUG_ONLY( std::cout << "dipole: polarization-angle: " << polarization_angle << std::endl; )
	DEBUG_ONLY( std::cout << "dipole: e-field-base: " << e_field_base << std::endl; )
	DEBUG_ONLY( std::cout << "dipole: e-field: " << e_field << std::endl; )

	// compute dipole vector (µS)
	for (int_t i = 0; i < sites; ++i)
		dipole_vector[i] = scalar_product(dipole_directions[i], e_field);

	return dipole_vector;
}

complex_matrix_t compute_dipole_matrix(const real_vector_t& dipole_vector,
                                       sites_to_states_mode_t mode)
{
	if (mode == sites_to_states_mode_t::identity)
		throw std::runtime_error("heom::compute_dipole_matrix(): sites_to_states_mode_t cannot be identity here due to physics.");

	const int_t sites = dipole_vector.size();
	const int_t states = sites_to_states(sites, mode);

	complex_matrix_t dipole_matrix(states, states, 0.0);
	// fill dipole_matrix matrix' first row from column 1 to sites (column 0 is ground state)
	for (int_t i = 0; i < sites; ++i)
		dipole_matrix.at(0, i + 1) = dipole_vector[i];

	// fill rest of dipole_matrix matrix
	if (mode == sites_to_states_mode_t::excited_state_absorption) {
		int_t n = sites + 1;
		for (int_t i = 0; i < sites; ++i)
			for (int_t j = i + 1; j < sites; ++j) {
				dipole_matrix.at(i + 1, n) = dipole_vector[j];
				dipole_matrix.at(j + 1, n) = dipole_vector[i];
				n++;
			}
	}

	return dipole_matrix;
}

} // namespace heom