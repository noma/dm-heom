// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_circular_dichroism_config_hpp
#define heom_circular_dichroism_config_hpp

#include <cstdint>
#include <fstream>
#include <boost/program_options.hpp>

#include "heom/common.hpp"
#include "heom/dipole_config_entries.hpp"
#include "heom/laser_config_entries.hpp"
#include "heom/config.hpp"

namespace heom {

/**
 * This class contains all options needed by the circular dichroism application
 */
class circular_dichroism_config : public config, public dipole_config_entries
{
public:
	/** Construct from file.
	 * config_file_name .. path to config file
	 * parse .. set to false to avoid automatic parsing (i.e. when called by a subclass)
	 */
	circular_dichroism_config(const std::string& config_file_name, bool parse_config = true);

	virtual ~circular_dichroism_config(){};

	// config value getters
	// NOTE: none

protected:
	virtual void check();

	virtual void check_final();

	virtual void post_process();

private:
	// config values
	// NOTE: none

};

} // namespace heom

#endif //  heom_circular_dichroism_config_hpp
