// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_static_fluorescence_config_hpp
#define heom_static_fluorescence_config_hpp

#include <cstdint>
#include <fstream>
#include <boost/program_options.hpp>

#include "heom/common.hpp"
#include "heom/dipole_config_entries.hpp"
#include "heom/laser_config_entries.hpp"
#include "heom/thermal_state_search_config.hpp"

namespace heom {

/**
 * This class contains all options needed by the fluorescence application
 */
class static_fluorescence_config : public thermal_state_search_config, public dipole_config_entries
{
public:
	/** Construct from file.
	 * config_file_name .. path to config file
	 * parse .. set to false to avoid automatic parsing (i.e. when called by a subclass)
	 */
	static_fluorescence_config(const std::string& config_file_name, bool parse_config = true);

	virtual ~static_fluorescence_config() {};

	// config value getters

protected:
	virtual void check();

	virtual void check_final();

	virtual void post_process();

private:
	// config values

};

} // namespace heom

#endif // heom_static_fluorescence_config_hpp
