// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_config_hpp
#define heom_config_hpp

#include <array>
#include <cstdint>
#include <fstream>

#include <noma/num/stepper_type.hpp>
#include <noma/typa/typa.hpp>

#include "heom/common.hpp"
#include "heom/config_base.hpp"
#include "heom/config_error.hpp"
#include "heom/filtering_strategy.hpp"
#include "heom/observation_type.hpp"

namespace heom {

/**
 * This class contains all options needed by all HEOM-applications.
 * Special needs are addressed by sub-classes that add further options.
 */
class config : public config_base
{
public:
	using obervation_type_list = parser::vector_wrapper<parser::pair_wrapper<observation_type_t, std::string>>;

	/** Construct from file.
	 * config_file_name .. path to config file
	 * parse .. set to false to avoid automatic parsing (i.e. when called by a subclass)
	 */
	config(const std::string& config_file_name, bool parse_config = true);

	virtual ~config() {};

	// config value getters

	// program
	const obervation_type_list& observations() const
	{ return observations_; }

	const int_t& program_observe_steps() const
	{ return program_observe_steps_; }

	// filter
	const heom::filtering_strategy_t& filtering_strategy() const
	 { return filtering_strategy_;}

	const int_t& filtering_first_layer() const
	 { return filtering_first_layer_; }

	// solver
	const num::stepper_type_t& solver_stepper_type() const
	{ return solver_stepper_type_; }

	const real_t& solver_step_size() const
	{ return solver_step_size_; }

	const int_t& solver_steps() const
	{ return solver_steps_; }

	const bool& solver_track_flows() const
	{ return solver_track_flows_; }

	const std::string& solver_flow_filename() const
	{ return solver_flow_filename_; }

	// system
	const int_t& system_sites() const
	{ return system_sites_; }

	const int_t& system_ado_depth() const
	{ return system_ado_depth_; }

	const complex_matrix_t& system_hamiltonian() const
	{ return system_hamiltonian_; }

	// baths
	const int_t& baths_number() const
	{ return baths_number_; }

	const int_t& baths_max_per_site() const
	{ return baths_max_per_site_; }

	const parser::matrix<int_t>& baths_coupling() const
	{ return baths_coupling_; }

	const real_vector_t& baths_lambda() const
	{ return baths_lambda_; }

	const real_vector_t& baths_invnu() const
	{ return baths_invnu_; }

	const real_vector_t& baths_uppercase_omega() const
	{ return baths_uppercase_omega_; }

	const int_t& baths_matsubaras() const
	{ return baths_matsubaras_; }

	const double& baths_temperature() const
	{ return baths_temperature_; }

protected:
	virtual void check();

	virtual void check_final();

	virtual void post_process();

private:
	// config values

	// program
	obervation_type_list observations_;
	int_t program_observe_steps_;

	// filtering
	// TODO: add two new entries: strategy, first_layer
	// if we choose filter strategy none we still have to set FILTER_FIRST_LAYER ..... ??!!
	// need some wrapper thing here?
	heom::filtering_strategy_t filtering_strategy_;
	int_t filtering_first_layer_;

	// solver
	num::stepper_type_t solver_stepper_type_;
	real_t solver_step_size_;
	int_t solver_steps_;
	bool solver_track_flows_;
	std::string solver_flow_filename_;

	// system
	int_t system_sites_;
	int_t system_ado_depth_;
	complex_matrix_t system_hamiltonian_;

	// baths
	int_t baths_number_;
	int_t baths_max_per_site_;
	parser::matrix<int_t> baths_coupling_;
	real_vector_t baths_lambda_;
	real_vector_t baths_invnu_;
	real_vector_t baths_uppercase_omega_;
	int_t baths_matsubaras_;
	real_t baths_temperature_;
};

} // namespace heom

#endif // heom_config_hpp
