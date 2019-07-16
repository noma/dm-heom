// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/static_fluorescence_config.hpp"

namespace heom {

namespace bpo = ::boost::program_options;

static_fluorescence_config::static_fluorescence_config(const std::string& config_file_name, bool parse_config)
	: thermal_state_search_config(config_file_name, false),
	  dipole_config_entries(desc_)
{
	// NOTE: no additional entries needed

	if (parse_config)
		parse(config_file_name);
}

void static_fluorescence_config::check()
{
	thermal_state_search_config::check(); // call direct super-class' check()

	// call generic checks from entry classes
	dipole_config_entries::check(system_sites());
}

void static_fluorescence_config::check_final()
{
	// check program_task_
	if (program_task() != program_task_t::STATIC_FLUORESCENCE)
		throw config_error("wrong program_task for this application");

	// check observations specification, static_fluorescence is a special case
	// as it averages multiple observations
	if (observations().get().size() > 1)
		throw config_error("static_fluorescence only supports a single observation of type matrix_trace_static_fluorescence");
	for (auto& obs : observations().get())
	{
		if (obs.get().first != heom::observation_type_t::matrix_trace_static_fluorescence)
			throw config_error("invalid observation type found in observations, only matrix_trace_static_fluorescence is allowed");
	}
}

void static_fluorescence_config::post_process()
{
	thermal_state_search_config::post_process(); // call direct super-class' post_process()

	// call post_process from entry classes
	dipole_config_entries::post_process();
}

} // namespace heom
