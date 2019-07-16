// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_laser_config_entries_hpp
#define heom_laser_config_entries_hpp

#include <boost/program_options.hpp>
#include <noma/typa/typa.hpp>

#include "heom/common.hpp"

namespace heom {

/**
 * This class contains the dipole entries, and can be used as super-class by
 * different app config classes to add dipole entries.
 */
class laser_config_entries
{
public:
	using laser_direction_type = real_vector_t;
	using laser_directions_type = parser::vector<laser_direction_type>;

	/** Construct from file.
	 * config_file_name .. path to config file
	 * parse .. set to false to avoid automatic parsing (i.e. when called by a subclass)
	 */
	laser_config_entries(boost::program_options::options_description& desc);

	// config value getters
	const laser_directions_type& laser_directions() const
	{ return laser_directions_; }

protected:
	void check();

	void post_process();

protected:
	// config values
	laser_directions_type laser_directions_;
};

} // namespace heom

#endif // heom_laser_config_entries_hpp
