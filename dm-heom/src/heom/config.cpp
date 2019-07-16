// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/config.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/hermitian.hpp>

#include "heom/constants.hpp"

namespace heom {

namespace bpo = ::boost::program_options;
namespace blas = ::boost::numeric::ublas;

config::config(const std::string& config_file_name, bool parse_config)
	: config_base(config_file_name, false)
{
	std::stringstream observation_type_help;
	observation_type_help << "List of observations, pairs {(observation_type, filename), ...}. Valid observation_types are: ";
	parser::write_map_value_list(observation_type_names, observation_type_help);

	// generate help text for filtering strategy
	std::stringstream filtering_strategy_help;
	filtering_strategy_help << "Strategy by which the hierarchy graph is filtered, one of: ";
	parser::write_map_value_list(filtering_strategy_names, filtering_strategy_help);

	std::stringstream solver_stepper_type_help;
	solver_stepper_type_help << "Stepper type for numerical integration, one of: ";
	parser::write_map_value_list(num::stepper_type_names, solver_stepper_type_help);

	// initialise program options description
	desc_.add_options()
		// program
		("program.observations", bpo::value(&observations_)->required(), observation_type_help.str().c_str())
		("program.observe_steps", bpo::value(&program_observe_steps_)->required(), "Observe (i.e. generate ouput) every n steps.")

		// filter
		("filtering.strategy", bpo::value(&filtering_strategy_)->required(), filtering_strategy_help.str().c_str())
		("filtering.first_layer", bpo::value(&filtering_first_layer_)->required(), "Hierarchy layer where filtering begins.")

		// solver
		("solver.stepper_type", bpo::value(&solver_stepper_type_)->required(), solver_stepper_type_help.str().c_str())
		("solver.step_size", bpo::value(&solver_step_size_)->required(), "Simulation time step in [s]")
		("solver.steps", bpo::value(&solver_steps_)->required(), "Simulated time [steps].")
		("solver.track_flows", bpo::value(&solver_track_flows_)->required(), "Enable flow tracking [true or false].")
		("solver.flow_filename", bpo::value(&solver_flow_filename_)->required(), "Output file for flow tracking.")

		// system
		("system.sites", bpo::value(&system_sites_)->required(), "Number of sites in the input system.")
		("system.ado_depth", bpo::value(&system_ado_depth_)->required(), "Depth of the hierarchy in the input system.")
		("system.hamiltonian", bpo::value(&system_hamiltonian_)->required(), "Hamiltonian, complex matrix (num-sites x num-sites).")

		// baths 
		("baths.number", bpo::value(&baths_number_)->required(), "Total number of independent baths.")
		("baths.max_per_site", bpo::value(&baths_max_per_site_)->required(), "Maximumexception number of baths coupled to a single site.")
		("baths.coupling", bpo::value(&baths_coupling_)->required(), "Baths connections to sites (baths_max_persite x system_sites), -1 for n/a.")
		("baths.lambda", bpo::value(&baths_lambda_)->required(), "Bath reorganization energies [lambda] for each bath [invcm].")
		("baths.invnu", bpo::value(&baths_invnu_)->required(), "Bath correlation times [nu] for each bath [inv fs].")
		("baths.Omega", bpo::value(&baths_uppercase_omega_)->required(), "Bath peak shift [Omega] for each bath [invcm].")
		("baths.matsubaras", bpo::value(&baths_matsubaras_)->required(), "Number of Matsubara frequencies per bath.")
		("baths.temperature", bpo::value(&baths_temperature_)->required(), "Temperature [K].")
		;

	if (parse_config)
		parse(config_file_name);
}

void config::check()
{
	config_base::check(); // call direct super-class' check()

	// observations
	if (observations().get().empty())
		throw config_error("observations must at least contain one entry.");

	DEBUG_ONLY( std::cout << "heom::config_base::parse(): observations=" << observations_ << std::endl; )

	// program_observe_steps
	if (!(program_observe_steps_ > 0))
		throw config_error("program.observe_steps must be > 0");

	// filtering_strategy
	if (filtering_strategy_ == filtering_strategy_t::NONE && filtering_first_layer_ != -1){
		throw config_error("filtering.first_layer must be -1, because it has no effect for filtering_strategy = none");
	}

	// filtering_first_layer
	if (filtering_strategy_ == filtering_strategy_t::SINGLE_EXCITATION && !(filtering_first_layer_ >= 0))
		throw config_error("filtering.first_layer must be greater or equal to zero for filtering.strategy = single excitation filter");

	// compare ADO cut off depth with filtering first layer
	if(system_ado_depth_ < filtering_first_layer_){
		throw config_error("filtering.first_layer cannot be larger than system_ado_depth");
	}

	// solver_step_size
	if (!(solver_step_size_ > 0.0)) // TODO: upper bound ?
		throw config_error("solver.step_size must be > 0");

	// solver_steps
	if (!(solver_steps_ > 0))
		throw config_error("solver.steps must be > 0");

	auto remainder_steps = solver_steps_ % program_observe_steps_;
	if (remainder_steps != 0)
		throw config_error("there's a remainder when dividing solver.steps by program.observe_steps");

	// solver_flow_filename
	// if (solver_track_flows_)
	// TODO: what happens if the path does not exist, what happens if the file exists already

	// system_sites
	if (!(system_sites_ > 0)) // TODO: upper bound
		throw config_error("system.sites must be > 0");

	// system_ado_depth
	if (!(system_ado_depth_ > 0)) // TODO: upper bound
		throw config_error("system.ado_depth must be > 0");

	// system_hamiltonian
	if (!(static_cast<int_t>(system_hamiltonian_.rows()) == system_sites_ &&
	      static_cast<int_t>(system_hamiltonian_.cols()) == system_sites_))
		throw config_error("system.hamiltonien must be a system.sites by system.sites matrix");

	// create boost ublas matrix from hamiltonian
	using matrix_t = blas::matrix<complex_t>;
	matrix_t hamiltonian_m(system_sites_, system_sites_);
	for (int_t i = 0; i < system_sites_; ++i)
		for (int_t j = 0; j < system_sites_; ++j)
			hamiltonian_m(i,j) = system_hamiltonian_.at(i, j);

	if(!blas::is_hermitian(hamiltonian_m))
		throw config_error("system.hamiltonien must be hermitian");

	// baths_number
	if (!(baths_number_ > 0)) // TODO: upper bound
		throw config_error("baths.number must be > 0");

	// baths_max_per_site
	if (!(baths_max_per_site_ > 0)) // TODO: upper bound
		throw config_error("baths.max_per_site must be > 0");
	if (!(baths_max_per_site_ <= baths_number_)) // TODO: upper bound
		throw config_error("baths.max_per_site must be <= baths.number");

	// baths_coupling
	if (!(static_cast<int_t>(baths_coupling_.rows()) == system_sites_ &&
	      static_cast<int_t>(baths_coupling_.cols()) == baths_max_per_site_))
		throw config_error("baths.coupling must be a system.sites by baths_max_per_site_ matrix");

	for (size_t i = 0; i < baths_coupling_.rows(); ++i)
		for (size_t j = 0; j < baths_coupling_.cols(); ++j)
			if (!(baths_coupling_.at(i,j) < baths_number_))
				throw config_error("baths.coupling values must between 0 and baths.number, or negative for non-coupling");

	// baths_lambda
	if (!(static_cast<int_t>(baths_lambda_.size()) == baths_number_))
		throw config_error("baths.lambda must have a length of baths.number");

	// baths_invnu
	if (!(static_cast<int_t>(baths_invnu_.size()) == baths_number_))
		throw config_error("baths.invnu must have a length of baths.number");
	for (size_t i = 0; i < baths_invnu_.size(); ++i)
		if (!(baths_invnu_.at(i) > 0)) // TODO: upper bound
			throw config_error("baths.invnu must have values > 0");

	// baths_uppercase_omage
	if (!(static_cast<int_t>(baths_uppercase_omega_.size()) == baths_number_))
		throw config_error("baths.Omega must have a length of baths.number");
	// TODO: if one Omega != 0, there must be a counterpart with same invnu, lambda, -Omega
}

void config::check_final()
{

}

void config::post_process()
{
	config_base::post_process(); // call direct super-class' post_process()

	// scale to SI units
	system_hamiltonian_.scale(constants::invcm_to_joule);
	baths_lambda_.scale(constants::invcm_to_joule);
	baths_uppercase_omega_.scale(constants::invcm_to_joule/constants::h_bar);
	baths_invnu_.scale(1.0E-15); // TODO: maybe change input format to SI seconds
}

} // namespace heom
