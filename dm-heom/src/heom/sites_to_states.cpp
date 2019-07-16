// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/sites_to_states.hpp"

#include <cmath>
#include <utility>
#include <vector>

namespace heom {

int_t sites_to_states(int_t sites, sites_to_states_mode_t mode)
{
	switch (mode) {
		case sites_to_states_mode_t::identity:
			return sites;
		case sites_to_states_mode_t::with_ground_state:
			return 1 + sites;
		case sites_to_states_mode_t::excited_state_absorption:
			return 1 + sites + ((sites * (sites - 1)) / 2); // NOTE: integeger division, but mathematically guaranteed to be fully divisible as on factor will be even
		default:
			throw std::runtime_error("heom::sites_to_states(): unhandled sites_to_states_mode_t value detected.");
	}
}

int_t max_baths_per_state(int_t max_baths_per_site, sites_to_states_mode_t mode)
{
	switch (mode) {
		case sites_to_states_mode_t::identity:
			return max_baths_per_site;
		case sites_to_states_mode_t::with_ground_state:
			return max_baths_per_site;
		case sites_to_states_mode_t::excited_state_absorption:
			return 2 * max_baths_per_site; // NOTE: lazy upper bound, actual value could be smaller
		default:
			throw std::runtime_error("heom::sites_to_states(): unhandled sites_to_states_mode_t value detected.");
	}
}

complex_matrix_t state_hamiltonian_with_ground_state(const complex_matrix_t& hamiltonian, int_t states)
{
	// zero initialised result matrix
	complex_matrix_t result(states, states, 0.0);

	// generate state hamiltonian
	// first row/columan stay zero, copy hamiltonian over to the remainder
	for (size_t i = 0; i < hamiltonian.rows(); ++i)
		for (size_t j = 0; j < hamiltonian.cols(); ++j)
			result.at(i + 1, j + 1) = hamiltonian.at(i, j);

	return result;
}

complex_matrix_t state_hamiltonian_excited_state_absorption(const complex_matrix_t& hamiltonian, int_t states)
{
	// zero initialised result matrix
	const int_t sites = hamiltonian.rows();
	complex_matrix_t result(states, states, 0.0);

	// generate state hamiltonian

	// copy over original hamiltonian next to the ground state row/col
	for (size_t i = 0; i < hamiltonian.rows(); ++i)
		for (size_t j = 0; j < hamiltonian.cols(); ++j)
			result.at(i + 1, j + 1) = hamiltonian.at(i, j);

	// create look up table for each state and its corresponding sites
	using int_pair = std::pair<int_t, int_t>;
	std::vector<int_pair> lut;
	lut.reserve(states);

	lut.push_back(int_pair {-1, -1});
	for (int_t i = 0; i < sites; ++i)
		lut.push_back(int_pair {i, i});

	// add all unique combination of sites
	for (int_t i = 0; i < sites; ++i)
		for (int_t j = i + 1; j < sites; ++j) {
			lut.push_back(int_pair {i, j});
		}

	// return 1 if all inpute are equal, zero otherwise
	auto kronecker_delta = [](int_t i, int_t j) { return i == j ? 1 : 0; };

	// fill rest of state_hamiltonian using LUT
	for (int_t i = sites + 1; i < states; ++i) {
		int_t ii = lut.at(i).first;
		int_t jj = lut.at(i).second;

		for (int_t j = sites + 1; j < states; ++j) {
			int_t kk = lut.at(j).first;
			int_t ll = lut.at(j).second;
			// NOTE: formula/scheme from some paper
			result.at(i, j) = static_cast<real_t>(kronecker_delta(ii, kk) * (1 - kronecker_delta(jj, ll))) * hamiltonian.at(jj, ll)
			                + static_cast<real_t>(kronecker_delta(ii, ll) * (1 - kronecker_delta(jj, kk))) * hamiltonian.at(jj, kk)
			                + static_cast<real_t>(kronecker_delta(jj, kk) * (1 - kronecker_delta(ii, ll))) * hamiltonian.at(ii, ll)
			                + static_cast<real_t>(kronecker_delta(jj, ll) * (1 - kronecker_delta(ii, kk))) * hamiltonian.at(ii, kk);
			// fix negative zeros that somehow occur in real part of
			if (result.at(i, j).real() == 0.0 && std::signbit(result.at(i, j).real())) {
				DEBUG_ONLY( std::cout << "heom::state_hamiltonian_excited_state_absorption(): found negative zero at states hamiltonian row " << i << " col " << j << ", fixing." << std::endl );
				result.at(i, j) = complex_t { 0.0, result.at(i, j).imag() };
			}
		}
		// diagonal
		result.at(i, i) = hamiltonian.at(ii, ii) + hamiltonian.at(jj, jj);
	}

	return result;
}

complex_matrix_t state_hamiltonian(const complex_matrix_t& hamiltonian, sites_to_states_mode_t mode)
{
	const int_t states = sites_to_states(hamiltonian.rows(), mode);

	switch (mode) {
		case sites_to_states_mode_t::identity:
			return hamiltonian;
		case sites_to_states_mode_t::with_ground_state:
			return state_hamiltonian_with_ground_state(hamiltonian, states);
		case sites_to_states_mode_t::excited_state_absorption:
			return state_hamiltonian_excited_state_absorption(hamiltonian, states);
		default:
			throw std::runtime_error("heom::state_hamiltonian(): unhandled sites_to_states_mode_t value detected.");
	}
}

parser::matrix<int_t> state_baths_coupling_with_ground_state(const parser::matrix<int_t>& coupling, int_t states)
{
	// new coupling matrix with a prepended ground state, initialised to -1, i.e. no coupling
	const int_t sites = coupling.rows();
	parser::matrix <int_t> result(sites + 1, coupling.cols(), -1);

	// copy original coupling
	for (size_t i = 1; i < result.rows(); ++i) // skip first row, i.e. ground state
		for (size_t j = 0; j < coupling.cols(); ++j)
			result.at(i, j) = coupling.at(i - 1, j); // one shift in rows

	return result;
}

parser::matrix<int_t> state_baths_coupling_excited_state_absorption(const parser::matrix<int_t>& coupling, int_t states)
{
	// new coupling matrix with a prepended ground state, and post-pended combined states
	const int_t sites = coupling.rows();
	const int_t with_ground_state_states = sites_to_states(sites, sites_to_states_mode_t::with_ground_state);

	parser::matrix<int_t> result(states, max_baths_per_state(coupling.cols(), sites_to_states_mode_t::excited_state_absorption), -1);

	// copy original coupling
	for (int_t i = 1; i < with_ground_state_states; ++i) // skip first row, i.e. ground state
		for (size_t j = 0; j < coupling.cols(); ++j)
			result.at(i, j) = coupling.at(i - 1, j); // one shift in rows

	// add combined coupling
	int_t k = with_ground_state_states;
	for (int_t i = 0; i < sites; ++i)
		for (int_t j = i + 1; j < sites; ++j) {
			// fill k-th row with union of baths coupled with state i and j

			// copy all couplings of site i
			size_t r = 0; // next result index to write
			for (size_t ib = 0; ib < coupling.cols(); ++ib) {
				if (coupling.at(i, ib) != -1) { // ignore -1
					result.at(k, r) = coupling.at(i, ib);
					++r;
				}
			}

			// copy all new couplings (we build a set union) of site j
			for (size_t jb = 0; jb < coupling.cols(); ++jb) {
				// check if current value is already present in result row k
				bool found = false;
				for (size_t rr = 0; rr < r; ++rr) { // iterate over alread written result
					if (coupling.at(j, jb) != -1) { // ignore -1
						if (result.at(k, rr) == coupling.at(j, jb)) {
							found = true;
							break;
						}

						if (!found) {
							result.at(k, r) = coupling.at(j, jb);
							++r;
						}
					}
				}
			}

			++k; // increment row
		} // for j

	return result;
}

parser::matrix<int_t> state_baths_coupling(const parser::matrix<int_t>& coupling, sites_to_states_mode_t mode)
{
	// rows is the number of sites, cols max_baths_per_site
	const int_t states = sites_to_states(coupling.rows(), mode);

	switch (mode) {
		case sites_to_states_mode_t::identity:
			return coupling;
		case sites_to_states_mode_t::with_ground_state:
			return state_baths_coupling_with_ground_state(coupling, states);
		case sites_to_states_mode_t::excited_state_absorption:
			return state_baths_coupling_excited_state_absorption(coupling, states);
		default:
			throw std::runtime_error("heom::state_baths_coupling(): unhandled sites_to_states_mode_t value detected.");
	}
}

} // namespace heom
