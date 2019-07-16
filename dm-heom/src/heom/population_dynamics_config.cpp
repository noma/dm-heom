// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/population_dynamics_config.hpp"

#include <fstream>
#include <iostream>

#include "heom/constants.hpp"

namespace heom {

namespace bpo = ::boost::program_options;

population_dynamics_config::population_dynamics_config(const std::string& config_file_name, bool parse_config)
	: config(config_file_name, false)

{
	// add options to boost program options description
	desc_.add_options()
		("population_dynamics.rho_init", bpo::value(&population_dynamics_rho_init_)->required(), "Initial density matrix (ONLY for task population_dynamics).")
		;

	if (parse_config)
		parse(config_file_name);
}

void population_dynamics_config::check()
{
	config::check(); // call direct super-class' check()

	// population_dynamics_rho_init
	if (!(static_cast<int_t>(population_dynamics_rho_init_.rows()) == system_sites() &&
	      static_cast<int_t>(population_dynamics_rho_init_.cols()) == system_sites()))
		throw config_error("population_dynamics.rho_init must be a system.sites by system.sites matrix");

}

void population_dynamics_config::check_final()
{
	// check program_task_
	if (program_task() != program_task_t::POPULATION_DYNAMICS)
		throw config_error("wrong program_task for this application");
}

void population_dynamics_config::post_process()
{
	config::post_process(); // call direct super-class' post_process()
}

} // namespace heom
