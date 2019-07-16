// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
// Copyright (c) 2016-2017 Lucas Deecke, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/thermal_state_search_config.hpp"

namespace heom {

namespace bpo = ::boost::program_options;

thermal_state_search_config::thermal_state_search_config(const std::string& config_file_name, bool parse_config)
	: config(config_file_name, false)
{
	// add options to boost program options description
	desc_.add_options()
		("thermal_state_search.delta", bpo::value(&thermal_state_search_delta_)->required(), "Convergence criteria: difference between hierarchy top norms between steps.")
		("thermal_state_search.max_steps", bpo::value(&thermal_state_search_max_steps_)->required(), "Abort criteria: stop solver after that many steps, set to 0 for unlimited.")
		;

	if (parse_config)
		parse(config_file_name);
}

void thermal_state_search_config::check()
{
	config::check(); // call direct super-class' check()
}

void thermal_state_search_config::check_final()
{
	// check program_task_
	if (program_task() != program_task_t::THERMAL_STATE_SEARCH)
		throw config_error("wrong program_task for this application");

	if (solver_step_size() != 1.0)
		throw config_error("solver.step_size_ must equal 1.0");

}

void thermal_state_search_config::post_process()
{
	config::post_process(); // call direct super-class' post_process()
}

} // namespace heom
