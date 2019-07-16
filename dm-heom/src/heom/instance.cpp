// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/instance.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include <boost/functional/hash.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include <noma/memory/memory.hpp>

#include "heom/constants.hpp"

namespace heom {

// constructor from heom_config
instance::instance(const config& cfg, sites_to_states_mode_t sts_mode, const hierarchy_graph& graph)
  : sites_to_states_mode_(sts_mode)
{
	// setup sizes
	states_ = sites_to_states(cfg.system_sites(), sts_mode); // NOTE: site to state conversion
	matsubaras_ = cfg.baths_matsubaras();
	baths_per_state_ = max_baths_per_state(cfg.baths_max_per_site(), sts_mode); // NOTE: site to state conversion
	baths_ = cfg.baths_number();

	ado_width_ = graph.width(); //baths_ * matsubaras_;
	ado_depth_ = graph.depth(); //cfg.system_ado_depth();

	matrices_ = graph.nodes(); // includes halo nodes
	first_non_plus_id_ = graph.first_non_plus_uid();
	assert(first_non_plus_id_ > 0);

	// allocate buffer memory
	allocate();

	// initialise low-level instance data structures from high_level graph partition data structure
	for (uid_t id = 0; id < graph.nodes(); ++id) {
		// copy tuples
		auto tuples = graph.tuples()[id].data();
		std::copy(tuples, tuples + ado_width_, &ado_tuples_[id * ado_width_]);
		// copy plus edges
		auto plus_edges = graph.plus_edges()[id].data();
		std::copy(plus_edges, plus_edges + ado_width_, &ado_plus_[id * ado_width_]);
		// copy minus edges
		auto minus_edges = graph.minus_edges()[id].data();
		std::copy(minus_edges, minus_edges + ado_width_, &ado_minus_[id * ado_width_]);
	}

	// initialise memory
	// NOTE: site to state conversion for the hamiltonian happens here
	set_hamiltonian(reinterpret_cast<const real_t*>(state_hamiltonian(cfg.system_hamiltonian(), sites_to_states_mode_).data()));

	// compute kernel input
	DEBUG_ONLY( std::cout << "instance::instance(..): begin precompute_baths_input" << std::endl; )
	// TODO: this call is annoyingly slow when testing
	precompute_baths_input(cfg);
	DEBUG_ONLY( std::cout << "instance::instance(..): end precompute_baths_input" << std::endl; )

	// initialize y_n (i.e. the initial value of our initial value problem)
	null_hierarchy();

	// NOTE: set_hierarchy_top() needs to be called to produce actual results, for testing: set (global) first element to 1.0
}

// constructor from saved state
instance::instance(const std::string& filename)
{
	allocate();
	read_file(filename);
}

// copy ctor
instance::instance(const instance& other)
{
	// copy shallow members
	version_ = other.version_;
	hierarchy_top_only_ = other.hierarchy_top_only_;
	hierarchy_layout_ = other.hierarchy_layout_;
	hamiltonian_layout_ = other.hamiltonian_layout_;
	sites_to_states_mode_ = other.sites_to_states_mode_;

	time_ = other.time_;
	first_non_plus_id_ = other.first_non_plus_id_;

	states_ = other.states_;
	matrices_ = other.matrices_;
	matsubaras_ = other.matsubaras_;
	baths_per_state_ = other.baths_per_state_;
	baths_ = other.baths_;
	ado_width_ = other.ado_width_;
	ado_depth_ = other.ado_depth_;

	allocate();

	// copy dynamic members
	std::copy(other.hierarchy_, other.hierarchy_ + other.size_hierarchy(), hierarchy_);
	std::copy(other.hamiltonian_, other.hamiltonian_ + other.size_hamiltonian(), hamiltonian_);

	std::copy(other.ado_tuples_, other.ado_tuples_ + other.size_ado_tuples(), ado_tuples_);
	std::copy(other.ado_plus_, other.ado_plus_ + other.size_ado_plus(), ado_plus_);
	std::copy(other.ado_minus_, other.ado_minus_ + other.size_ado_minus(), ado_minus_);

	std::copy(other.pf_same_sum_, other.pf_same_sum_ + other.size_pf_same_sum(), pf_same_sum_);
	std::copy(other.pf_same_, other.pf_same_ + other.size_pf_same(), pf_same_);
	std::copy(other.pf_circ_, other.pf_circ_ + other.size_pf_circ(), pf_circ_);
	std::copy(other.cbk_, other.cbk_ + other.size_cbk(), cbk_);
	std::copy(other.cbk_sqrt_, other.cbk_sqrt_ + other.size_cbk_sqrt(), cbk_sqrt_);
}

instance::~instance()
{
	free();
}

void instance::allocate()
{
	hierarchy_ = memory::aligned::instance().allocate<real_t>(size_hierarchy());
	hamiltonian_ = memory::aligned::instance().allocate<real_t>(size_hamiltonian());

	ado_tuples_ = memory::aligned::instance().allocate<std::int32_t>(size_ado_tuples());
	ado_plus_ = memory::aligned::instance().allocate<std::int32_t>(size_ado_plus());
	ado_minus_ = memory::aligned::instance().allocate<std::int32_t>(size_ado_minus());

	pf_same_sum_ = memory::aligned::instance().allocate<real_t>(size_pf_same_sum());
	pf_same_ = memory::aligned::instance().allocate<real_t>(size_pf_same());
	pf_circ_ = memory::aligned::instance().allocate<real_t>(size_pf_circ());
	cbk_ = memory::aligned::instance().allocate<real_t>(size_cbk());
	cbk_sqrt_ = memory::aligned::instance().allocate<real_t>(size_cbk_sqrt());
}

void instance::free()
{
	memory::aligned::instance().free(hierarchy_);
	memory::aligned::instance().free(hamiltonian_);

	memory::aligned::instance().free(ado_tuples_);
	memory::aligned::instance().free(ado_plus_);
	memory::aligned::instance().free(ado_minus_);

	memory::aligned::instance().free(pf_same_sum_);
	memory::aligned::instance().free(pf_same_);
	memory::aligned::instance().free(pf_circ_);
	memory::aligned::instance().free(cbk_);
	memory::aligned::instance().free(cbk_sqrt_);
}

void instance::precompute_baths_input(const config& cfg)
{

	complex_t* PFSAMESUMIn = reinterpret_cast<complex_t*>(this->pf_same_sum_);
	complex_t* PFSAMEIn = reinterpret_cast<complex_t*>(this->pf_same_);
	complex_t* PFCIRCIn = reinterpret_cast<complex_t*>(this->pf_circ_);
	complex_t* CBKIn = reinterpret_cast<complex_t*>(this->cbk_);
	real_t* CBKSqrtIn = reinterpret_cast<real_t*>(this->cbk_sqrt_);

	// FIXME: location of physical constants definitions?
	const real_t& hbar = constants::h_bar;

	real_t beta = 1.0 / (constants::k_b * cfg.baths_temperature());

	for(int_t b = 0; b < cfg.baths_number(); ++b)
	{
		complex_t aux;

		complex_t lambda_b { cfg.baths_lambda().at(b) }; // complexify
		complex_t nu_b { 1.0/cfg.baths_invnu().at(b), cfg.baths_uppercase_omega().at(b) }; // complexify

		PFCIRCIn[b] = nu_b * lambda_b / hbar;
		aux = 2.0 * lambda_b / (beta * hbar * hbar * nu_b);

		for (int_t k = 0; k < matsubaras_; ++k)
		{
			complex_t cbk, cbk0, nubk;
			// only the k=0 term differs with respect to complex conjugation
			if (k == 0)
			{
				//NOTE: with PADE
				complex_t auxsum=0.0;
				for(int_t j = 0; j < matsubaras_; ++j)
				{
					real_t eta = constants::eta_pade[(matsubaras_ - 1) * constants::max_pade + j];
					real_t xi = constants::xi_pade [(matsubaras_ - 1) * constants::max_pade + j];
					auxsum += (2.0 * eta * nu_b * nu_b) / (xi * xi / (beta * beta * hbar * hbar) - nu_b * nu_b);
				}
				cbk0 = 2.0 * lambda_b / (beta * hbar) * (1.0 - auxsum);

				// Alternative without PADE
				//cbk0  =lambda_b*nu_b*(1.0/tan(0.5*beta*hbar*nu_b));

				cbk = cbk0;
				nubk = nu_b;
			}
			else
			{

				//NOTE: with PADE
				real_t eta = constants::eta_pade[(matsubaras_ - 1) * constants::max_pade + k];
				real_t xi = constants::xi_pade[(matsubaras_ - 1) * constants::max_pade + k];
				cbk = 4.0 * lambda_b * nu_b / (beta * beta * hbar * hbar) * (eta * xi) / (xi * xi / (beta * beta * hbar * hbar) - nu_b * nu_b);
				nubk = xi / (beta*hbar);

				// Alternative without PADE
				//cbk=8.0*real_t(k)*M_PI*lambda_b*nu_b/(4.0*k*k*M_PI*M_PI-beta*beta*nu_b*nu_b*hbar*hbar);
				//nubk=2.0*real_t(k)*M_PI/(beta*hbar);

			}
			// FIXME we define a mapping here
			CBKIn[b * matsubaras_ + k] = cbk;
			CBKSqrtIn[b * matsubaras_ + k] = sqrt(abs(cbk));
			aux -= cbk / (hbar * nubk);
		}
		PFSAMEIn[b] = aux;
	}

	// pre-compute outer summation:
	// TODO: does not work with filtering, as it assumes all matrices are included in the propagation
	// TODO: do when coding post-filter compactify...

	for (uid_t id = 0; id < matrices_; ++id)
	{
		// long double needed due to large summation to avoid cancellation effects
		long double zr = 0.0L;
		long double zi = 0.0L;

		for (int b = 0; b < cfg.baths_number(); ++b)
		{
			complex_t nu_b { 1.0 / cfg.baths_invnu().at(b), cfg.baths_uppercase_omega().at(b) }; // complexify

			for (int k = 0; k < matsubaras_; ++k)
			{
				long double nubk_r,nubk_i;
				if (k==0)
				{
					nubk_r = static_cast<long double>(real(nu_b));
					nubk_i = static_cast<long double>(imag(nu_b));
				}
				else
				{
					// NOTE: with PADE
					nubk_r = constants::xi_pade[(matsubaras_ - 1) * constants::max_pade + k] / static_cast<long double>(beta*hbar);
					// Alternative without PADE
					//nubk_r=2.0L*(long double)M_PI*(long double)k/(beta*hbar);

					nubk_i = 0.0L;
				}
				// NOTE: we imply a mapping here
				int m = b * matsubaras_ + k;

				zr -= nubk_r * static_cast<long double>(ado_tuples_[id * ado_width_ + m]);
				zi -= nubk_i * static_cast<long double>(ado_tuples_[id * ado_width_ + m]);
			}
		}
		PFSAMESUMIn[id] = complex_t(zr, zi);
	} // sum over all HEOM matrices

	// output for comparison
//	DEBUG_ONLY( array_to_file("HAMIn", reinterpret_cast<complex_t*>(this->hamiltonian_), this->size_hamiltonian_byte() / 2 / sizeof(real_t)); )
//	DEBUG_ONLY( array_to_file("ADOTupleIn", this->ado_tuples_, this->size_ado_tuples_byte() / sizeof(std::int32_t)); )
//	DEBUG_ONLY( array_to_file("PlusIndexIn", this->ado_plus_, this->size_ado_plus_byte() / sizeof(std::int32_t)); )
//	DEBUG_ONLY( array_to_file("MinusIndexIn", this->ado_minus_, this->size_ado_minus_byte() / sizeof(std::int32_t)); )
//	DEBUG_ONLY( array_to_file("PFSAMESUMIn", PFSAMESUMIn, this->size_pf_same_sum_byte() / 2 / sizeof(real_t)); )
//	DEBUG_ONLY( array_to_file("PFSAMEIn", PFSAMEIn, this->size_pf_same_byte() / 2 / sizeof(real_t)); )
//	DEBUG_ONLY( array_to_file("PFCIRCIn", PFCIRCIn, this->size_pf_circ_byte() / 2 / sizeof(real_t)); )
//	DEBUG_ONLY( array_to_file("CBKIn", CBKIn, this->size_cbk_byte() / 2 / sizeof(real_t)); )
}

template<typename T>
std::ifstream::char_type* cast_to_char(T* ptr)
{
	return reinterpret_cast<std::ifstream::char_type*>(ptr);
}

bool instance::read_file(std::string filename)
{
	try {
		std::ifstream file;
		file.open(filename, std::ios::in | std::ios::binary);

		if (file.fail())
			throw config_error("could not open instance file: " + filename);

		// format specifics
		file.read(cast_to_char(&version_), sizeof version_);
		file.read(cast_to_char(&hierarchy_top_only_), sizeof hierarchy_top_only_);
		file.read(cast_to_char(&hierarchy_layout_), sizeof hierarchy_layout_);
		file.read(cast_to_char(&hamiltonian_layout_), sizeof hamiltonian_layout_);
		file.read(cast_to_char(&sites_to_states_mode_), sizeof sites_to_states_mode_);

		// variables
		file.read(cast_to_char(&time_), sizeof time_);
		file.read(cast_to_char(&first_non_plus_id_), sizeof first_non_plus_id_);

		// sizes
		file.read(cast_to_char(&states_), sizeof states_);
		file.read(cast_to_char(&matrices_), sizeof matrices_);
		file.read(cast_to_char(&matsubaras_), sizeof matsubaras_);
		file.read(cast_to_char(&baths_per_state_), sizeof baths_per_state_);
		file.read(cast_to_char(&baths_), sizeof baths_);
		file.read(cast_to_char(&ado_width_), sizeof ado_width_);
		file.read(cast_to_char(&ado_depth_), sizeof ado_depth_);

		// buffers
		file.read(cast_to_char(hierarchy_), hierarchy_top_only_ ? size_hierarchy_top_byte() : size_hierarchy_byte());
		file.read(cast_to_char(hamiltonian_), size_hamiltonian_byte());
		file.read(cast_to_char(ado_tuples_), size_ado_tuples_byte());
		file.read(cast_to_char(ado_plus_), size_ado_plus_byte());
		file.read(cast_to_char(ado_minus_), size_ado_minus_byte());
		file.read(cast_to_char(pf_same_sum_), size_pf_same_sum_byte());
		file.read(cast_to_char(pf_same_), size_pf_same_byte());
		file.read(cast_to_char(pf_circ_), size_pf_circ_byte());
		file.read(cast_to_char(cbk_), size_cbk_byte());
		file.read(cast_to_char(cbk_sqrt_), size_cbk_sqrt_byte());

		file.close();
	} catch (std::exception& e) {
		std::cerr << "instance::read_file(): Exception for file '" << filename << "': " << e.what() << std::endl;
		return false;
	}

	return true;

}

bool instance::write_file(std::string filename)
{
	try {
		std::ofstream file;
		file.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);

		// format specifics
		file.write(cast_to_char(&version_), sizeof version_);
		file.write(cast_to_char(&hierarchy_top_only_), sizeof hierarchy_top_only_);
		file.write(cast_to_char(&hierarchy_layout_), sizeof hierarchy_layout_);
		file.write(cast_to_char(&hamiltonian_layout_), sizeof hamiltonian_layout_);
		file.write(cast_to_char(&sites_to_states_mode_), sizeof sites_to_states_mode_);

		// variables
		file.write(cast_to_char(&time_), sizeof time_);
		file.write(cast_to_char(&first_non_plus_id_), sizeof first_non_plus_id_);

		// sizes
		file.write(cast_to_char(&states_), sizeof states_);
		file.write(cast_to_char(&matrices_), sizeof matrices_);
		file.write(cast_to_char(&matsubaras_), sizeof matsubaras_);
		file.write(cast_to_char(&baths_per_state_), sizeof baths_per_state_);
		file.write(cast_to_char(&baths_), sizeof baths_);
		file.write(cast_to_char(&ado_width_), sizeof ado_width_);
		file.write(cast_to_char(&ado_depth_), sizeof ado_depth_);

		// buffers
		file.write(cast_to_char(hierarchy_), hierarchy_top_only_ ? size_hierarchy_top_byte() : size_hierarchy_byte());
		file.write(cast_to_char(hamiltonian_), size_hamiltonian_byte());
		file.write(cast_to_char(ado_tuples_), size_ado_tuples_byte());
		file.write(cast_to_char(ado_plus_), size_ado_plus_byte());
		file.write(cast_to_char(ado_minus_), size_ado_minus_byte());
		file.write(cast_to_char(pf_same_sum_), size_pf_same_sum_byte());
		file.write(cast_to_char(pf_same_), size_pf_same_byte());
		file.write(cast_to_char(pf_circ_), size_pf_circ_byte());
		file.write(cast_to_char(cbk_), size_cbk_byte());
		file.write(cast_to_char(cbk_sqrt_), size_cbk_sqrt_byte());

		file.close();
	} catch (std::exception& e) {
		std::cerr << "instance::write_file(): Exception for file '" << filename << "': " << e.what() << std::endl;
		return false;
	}

	return true;
}

size_t instance::size_hierarchy() const
{
	return states_ * states_ * matrices_ * 2;
}

size_t instance::size_hierarchy_top() const
{
	return states_ * states_ * 2; // just the top of the hierarchy
}

size_t instance::size_hamiltonian() const
{
	return states_ * states_ * 2;
}

size_t instance::size_ado_tuples() const
{
	return matrices_ * ado_width_;
}

size_t instance::size_ado_plus() const
{
	return matrices_ * ado_width_;
}

size_t instance::size_ado_minus() const
{
	return matrices_ * ado_width_;
}

size_t instance::size_pf_same_sum() const
{
	return matrices_ * 2;
}

size_t instance::size_pf_same() const
{
	return baths_ * 2;
}

size_t instance::size_pf_circ() const
{
	return baths_ * 2;
}

size_t instance::size_cbk() const
{
	return baths_ * matsubaras_ * 2;
}

size_t instance::size_cbk_sqrt() const
{
	return baths_ * matsubaras_;
}

size_t instance::allocated_byte() const
{
	return size_hierarchy_byte() +
	       size_hamiltonian_byte() +
	       size_ado_tuples_byte() +
	       size_ado_plus_byte() +
	       size_ado_minus_byte() +
	       size_pf_same_sum_byte() +
	       size_pf_same_byte() +
	       size_pf_circ_byte() +
	       size_cbk_byte() +
	       size_cbk_sqrt_byte();
}

void instance::null_hierarchy() const
{
	for (size_t i = 0; i < this->size_hierarchy_byte() / sizeof(real_t); ++i) {
		hierarchy_[i] = 0.0;
	}
}

void instance::set_hamiltonian(const complex_matrix_t& new_hamiltonian) const
{
	set_hamiltonian(reinterpret_cast<const real_t*>(new_hamiltonian.data()));
}

void instance::set_hamiltonian(const real_t* new_hamiltonian) const
{
	for (size_t j = 0; j < size_hamiltonian(); ++j) {
		hamiltonian_[j] = new_hamiltonian[j];
	}
}

void instance::set_hamiltonian(const std::valarray<complex_t> new_hamiltonian) const
{
	for (size_t j = 0; j < size_hamiltonian(); j += 2) {
		hamiltonian_[j] = new_hamiltonian[j / 2].real(); // NOTE: intended integer division
		hamiltonian_[j+1] = new_hamiltonian[j / 2].imag();
	}
}

void instance::get_hierarchy_top(complex_t* buffer) const
{
	// TODO: adapt for when other memory layouts are introduced, maybe use some overloaded []
	for (int_t j = 0; j < states_ * states_; j += 1) {
		buffer[j] = complex_t( hierarchy_[j * 2], hierarchy_[j * 2 + 1] );
	}
}

void instance::set_hierarchy_top() const
{
	hierarchy_[0] = 1.0;
}

void instance::set_hierarchy_top(const complex_t* new_hierarchy_top) const
{
	// TODO: adapt for when other memory layouts are introduced, maybe use some overloaded []
	for (int_t j = 0; j < states_ * states_; ++j) {
		hierarchy_[2 * j] = new_hierarchy_top[j].real();
		hierarchy_[2 * j + 1] = new_hierarchy_top[j].imag();
	}
}

} // namespace heom
