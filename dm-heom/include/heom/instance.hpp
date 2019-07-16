// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_instance_hpp
#define heom_instance_hpp

#include <cstdint>
#include <string>
#include <valarray>

#include "noma/memory/array_to_file.hpp"

#include "heom/common.hpp"
#include "heom/config.hpp"
#include "heom/hierarchy_graph.hpp"
#include "heom/integer_partition.hpp"
#include "heom/sites_to_states.hpp"

namespace heom {

class hierarchy_graph;

class instance
{
public:
	// types
	enum class version_t : std::int32_t
	{
		V_1
	};
	typedef uint8_t byte_t;
	enum class memory_layout_t : std::int32_t
	{
		VANILLA
	};

	// version information
	static const version_t impl_version = version_t::V_1; // implementation version (this file)

	// constructor from heom_config
	instance(const config& config, sites_to_states_mode_t sts_mode, const hierarchy_graph& graph);

	// constructor from saved state
	instance(const std::string& filename);

	// copy ctor
	instance(const instance& other);

	~instance();

	instance& operator=(const instance& other) = delete;

	sites_to_states_mode_t sites_to_states_mode() const { return sites_to_states_mode_; }

	// i/o
	// TODO: maybe overload stream operators and use them here, for more general usage
	bool read_file(std::string filename); // returns false on error
	bool write_file(std::string filename); // returns false on error

	// buffer_size computations, results in sizeof(T)
	size_t size_hierarchy() const;
	size_t size_hierarchy_top() const;
	size_t size_hamiltonian() const;
	size_t size_ado_tuples() const;
	size_t size_ado_plus() const;
	size_t size_ado_minus() const;
	size_t size_pf_same_sum() const;
	size_t size_pf_same() const;
	size_t size_pf_circ() const;
	size_t size_cbk() const;
	size_t size_cbk_sqrt() const;

	// buffer size computations, results in byte
	size_t size_hierarchy_byte() const { return size_hierarchy() * sizeof(real_t); }
	size_t size_hierarchy_top_byte() const { return size_hierarchy_top() * sizeof(real_t); }
	size_t size_hamiltonian_byte() const { return size_hamiltonian() * sizeof(real_t); }
	size_t size_ado_tuples_byte() const { return size_ado_tuples() * sizeof(std::int32_t); }
	size_t size_ado_plus_byte() const { return size_ado_plus() * sizeof(std::int32_t); }
	size_t size_ado_minus_byte() const { return size_ado_minus() * sizeof(std::int32_t); }
	size_t size_pf_same_sum_byte() const { return size_pf_same_sum() * sizeof(real_t); }
	size_t size_pf_same_byte() const { return size_pf_same() * sizeof(real_t); }
	size_t size_pf_circ_byte() const { return size_pf_circ() * sizeof(real_t); }
	size_t size_cbk_byte() const { return size_cbk() * sizeof(real_t); }
	size_t size_cbk_sqrt_byte() const { return size_cbk_sqrt() * sizeof(real_t); }

	size_t allocated_byte() const;

	real_t time() const
	{ return time_; }

	std::int32_t first_non_plus_id() const
	{ return first_non_plus_id_; }

	std::int64_t states() const
	{ return states_; }

	std::int64_t matrices() const
	{ return matrices_; }

	std::int64_t matsubaras() const
	{ return matsubaras_; }

	std::int64_t baths_per_state() const
	{ return baths_per_state_; }

	std::int64_t baths() const
	{ return baths_; }

	std::int64_t ado_width() const
	{ return ado_width_; }

	std::int64_t ado_depth() const
	{ return ado_depth_; }

	real_t* hierarchy() const
	{ return hierarchy_; }

	const real_t* hamiltonian() const
	{ return hamiltonian_; }

	const std::int32_t* ado_tuples() const
	{ return ado_tuples_; }

	const std::int32_t* ado_plus() const
	{ return ado_plus_; }

	const std::int32_t* ado_minus() const
	{ return ado_minus_; }

	// change hierarchy after initialisation
	void null_hierarchy() const;

	// change hamiltonian after initialisation
	void set_hamiltonian(const complex_matrix_t& new_hamiltonian) const;
	void set_hamiltonian(const real_t* new_hamiltonian) const;
	void set_hamiltonian(const std::valarray<complex_t> new_hamiltonian) const;

	/**
	 * Get the top-level matrix in standard complex memory layout.
	 */
	void get_hierarchy_top(complex_t* buffer) const;

	/**
	 * Default initialise the top-level matrix by setting the first element to (1.0, 0.0).
	 */
	void set_hierarchy_top() const;

	/**
	 * Set the top-level matrix from standard complex memory layout.
	 */
	void set_hierarchy_top(const complex_t* new_hierarchy_top) const;


	const real_t* pf_same_sum() const
	{ return pf_same_sum_; }

	const real_t* pf_same() const
	{ return pf_same_; }

	const real_t* pf_circ() const
	{ return pf_circ_; }

	const real_t* cbk() const
	{ return cbk_; }

	const real_t* cbk_sqrt() const
	{ return cbk_sqrt_; }

private:
	// memory management 
	void allocate(); // all sizes should be set correctly!
	void free(); // all sizes should be set correctly!

	void precompute_baths_input(const config& cfg);

	// format specifics
	// TODO: all unused
	version_t version_ = impl_version; // file format version
	// TOOD: remove hierarchy_top_only?
	bool hierarchy_top_only_ = false; // only the first sigma matrix of the hierarchy is actually stored, everything else shall be null-initialised
	memory_layout_t hierarchy_layout_ = memory_layout_t::VANILLA;
	memory_layout_t hamiltonian_layout_ = memory_layout_t::VANILLA;
	sites_to_states_mode_t sites_to_states_mode_ = sites_to_states_mode_t::identity; // NOTE: some default

	// variables
	real_t time_ = 0.0; // TODO: unused
	std::int32_t first_non_plus_id_; // id of the first sigma-matrix without a plus link

	// sizes
	std::int64_t states_;
	std::int64_t matrices_;
	std::int64_t matsubaras_;
	std::int64_t baths_per_state_;
	std::int64_t baths_;
	std::int64_t ado_width_;
	std::int64_t ado_depth_;

	// buffers
	real_t* hierarchy_; // array of sigma matrices, complex: states_ * states_ * matrices_ * 2 * sizeof(real_t)
	real_t* hamiltonian_; // hamiltonian, complex: states_ * states_ * 2 * sizeof(real_t)

	std::int32_t* ado_tuples_; // ADO tuples: matrices_ * ado_width * sizeof(std::int32_t)
	std::int32_t* ado_plus_; // ADO plus indices: matrices_ * ado_width * sizeof(std::int32_t)
	std::int32_t* ado_minus_; // ADO minux indices: matrices_ * ado_width * sizeof(std::int32_t)

	real_t* pf_same_sum_; // complex: matrices_ * 2 * sizeof(real_t)
	real_t* pf_same_; // complex: bath * 2 * sizeof(real_t)
	real_t* pf_circ_; // complex: bath * 2 * sizeof(real_t)
	real_t* cbk_; // complex: bath * matsubaras_ * 2 * sizeof(real_t)
	real_t* cbk_sqrt_; // real: bath * matsubaras_ * sizeof(real_t)
}; // class instance

} // namespace heom

#endif // heom_instance_hpp
