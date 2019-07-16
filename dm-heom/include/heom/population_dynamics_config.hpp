// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_population_dynamics_config_hpp
#define heom_population_dynamics_config_hpp

#include <cstdint>
#include <fstream>

#include <boost/program_options.hpp>
#include <noma/typa/matrix.hpp>

#include "heom/common.hpp"
#include "heom/config.hpp"

namespace heom {

class population_dynamics_config : public config
{
public:
	/** Construct from file.
	 * config_file_name .. path to config file
	 * parse .. set to false to avoid automatic parsing (i.e. when called by a subclass)
	 */
	population_dynamics_config(const std::string& config_file_name, bool parse_config = true);

	virtual ~population_dynamics_config() {};

	// config value getters

	const complex_matrix_t& population_dynamics_rho_init() const
	{ return population_dynamics_rho_init_; }

protected:
	virtual void check();

	virtual void check_final();

	virtual void post_process();

private:
	// config values
	complex_matrix_t population_dynamics_rho_init_;
};

} // namespace heom

#endif // heom_population_dynamics_config_hpp
