// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_rate_kernel_config_hpp
#define heom_rate_kernel_config_hpp

#include <cstdint>
#include <fstream>

#include <boost/program_options.hpp>
#include <noma/typa/matrix.hpp>

#include "heom/common.hpp"
#include "heom/config.hpp"

namespace heom {

class rate_kernel_config : public config
{
public:
	/** Construct from file.
	 * config_file_name .. path to config file
	 * parse .. set to false to avoid automatic parsing (i.e. when called by a subclass)
	 */
	rate_kernel_config(const std::string& config_file_name, bool parse_config = true);

	virtual ~rate_kernel_config() {};

	// config value getters
	// NOTE: nothing needed here

protected:
	virtual void check();

	virtual void check_final();

	virtual void post_process();

private:
	// config values
	// NOTE: nothing needed here
};

} // namespace heom

#endif // heom_rate_kernel_config_hpp
