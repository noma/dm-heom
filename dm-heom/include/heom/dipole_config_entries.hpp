// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_dipole_config_entries_hpp
#define heom_dipole_config_entries_hpp

#include <boost/program_options.hpp>
#include <noma/typa/typa.hpp>

#include "heom/common.hpp"

namespace heom {

/**
 * This class contains the dipole entries, and can be used as super-class by
 * different app config classes to add dipole entries.
 */
class dipole_config_entries
{
public:
	using dipole_direction_type = real_vector_t;
	using dipole_directions_type = parser::vector<dipole_direction_type>;
	using dipole_strengths_type = parser::vector<real_t>;
	using dipole_centers_type = dipole_directions_type;
	using dipole_tensor_prefactors_type = real_vector_t;
	using dipole_tensor_components_type = parser::vector<parser::vector<int_t>>;

	/** Construct from file.
	 * config_file_name .. path to config file
	 * parse .. set to false to avoid automatic parsing (i.e. when called by a subclass)
	 */
	dipole_config_entries(boost::program_options::options_description& desc);

	// config value getters
	const dipole_directions_type& dipole_directions() const
	{ return dipole_directions_; }

	const dipole_strengths_type& dipole_strengths() const
	{ return dipole_strengths_; }

	const dipole_centers_type& dipole_centers() const
	{ return dipole_centers_; }

	const dipole_tensor_prefactors_type& dipole_tensor_prefactors() const
	{ return dipole_tensor_prefactors_; }

	const dipole_tensor_components_type& dipole_tensor_components() const
	{ return dipole_tensor_components_; }


protected:
	void check(int_t system_sites);

	void post_process();

protected:
	// config values
	dipole_directions_type dipole_directions_;
	dipole_strengths_type dipole_strengths_;
	dipole_centers_type dipole_centers_;
	dipole_tensor_prefactors_type dipole_tensor_prefactors_;
	dipole_tensor_components_type dipole_tensor_components_;
};

} // namespace heom

#endif // heom_dipole_config_entries_hpp
