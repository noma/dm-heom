// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
// Copyright (c) 2016-2017 Lucas Deecke, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_thermal_state_search_config_hpp
#define heom_thermal_state_search_config_hpp

#include <cstdint>
#include <fstream>

#include <boost/program_options.hpp>
#include <noma/typa/matrix.hpp>

#include "heom/common.hpp"
#include "heom/config.hpp"

namespace heom {

/**
 * This class contains all options needed by all HEOM-applications.
 * Special needs are addressed by sub-classes that add further options.
 */
class thermal_state_search_config : public config
{
public:
	/** Construct from file.
	 * config_file_name .. path to config file
	 * parse .. set to false to avoid automatic parsing (i.e. when called by a subclass)
	 */
	thermal_state_search_config(const std::string& config_file_name, bool parse_config = true);

	virtual ~thermal_state_search_config() {};

	// config value getters

	const real_t& thermal_state_search_delta() const
	{ return thermal_state_search_delta_; }

	const int_t& thermal_state_search_max_steps() const
	{ return thermal_state_search_max_steps_; }

protected:
	virtual void check();

	virtual void check_final();

	virtual void post_process();

private:
	// config values
	real_t thermal_state_search_delta_;
	int_t thermal_state_search_max_steps_;
};

} // namespace heom

#endif // heom_thermal_state_search_config_hpp
