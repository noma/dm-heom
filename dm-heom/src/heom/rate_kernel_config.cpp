// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/rate_kernel_config.hpp"

#include <fstream>
#include <iostream>

#include "heom/constants.hpp"

namespace heom {

namespace bpo = ::boost::program_options;

rate_kernel_config::rate_kernel_config(const std::string& config_file_name, bool parse_config)
	: config(config_file_name, false)

{
	// add options to boost program options description
	// NOTE: nothing needed here

	if (parse_config)
		parse(config_file_name);
}

void rate_kernel_config::check()
{
	config::check(); // call direct super-class' check()
}

void rate_kernel_config::check_final()
{
	// check program_task_
	if (program_task() != program_task_t::RATE_KERNEL)
		throw config_error("wrong program_task for this application");

	if (observations().get().size() > 1)
		throw config_error("rate_kernel only supports a single observation of type matrix_complete");
	for (auto& obs : observations().get())
	{
		std::cout << obs.get().first << std::endl;
		if (obs.get().first != heom::observation_type_t::matrix_complete_rate_kernel)
			throw config_error("invalid observation type found in observations, only matrix_complete_rate_kernel is allowed");
	}

}

void rate_kernel_config::post_process()
{
	config::post_process(); // call direct super-class' post_process()
}

} // namespace heom
