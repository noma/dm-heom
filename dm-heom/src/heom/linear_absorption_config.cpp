// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/linear_absorption_config.hpp"

#include <fstream>
#include <iostream>

#include "heom/constants.hpp"

namespace heom {

namespace bpo = ::boost::program_options;

linear_absorption_config::linear_absorption_config(const std::string& config_file_name, bool parse_config)
	: config(config_file_name, false),
	  dipole_config_entries(desc_)
{
	// add options to boost program options description
	// NOTE: none

	if (parse_config)
		parse(config_file_name);
}

void linear_absorption_config::check()
{
	config::check(); // call direct super-class' check()

	// call generic checks from entry classes
	dipole_config_entries::check(system_sites());
}

void linear_absorption_config::check_final()
{
	// check program_task_
	if (program_task() != program_task_t::LINEAR_ABSORPTION)
		throw config_error("wrong program_task for this application");
}

void linear_absorption_config::post_process()
{
	config::post_process(); // call direct super-class' post_process()

	// call post_process from entry classes
	dipole_config_entries::post_process();
}


} // namespace heom
