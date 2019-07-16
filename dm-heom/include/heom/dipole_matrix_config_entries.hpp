// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_dipole_matrix_config_entries_hpp
#define heom_dipole_matrix_config_entries_hpp

#include <boost/program_options.hpp>
#include <noma/typa/typa.hpp>

#include "heom/common.hpp"

namespace heom {

/**
 * This class contains the dipole entries, and can be used as super-class by
 * different app config classes to add dipole entries.
 */
class dipole_matrix_config_entries
{
public:
	using dipole_matrix_type = parser::vector<complex_matrix_t>;

	/** Construct from file.
	 * config_file_name .. path to config file
	 * parse .. set to false to avoid automatic parsing (i.e. when called by a subclass)
	 */
	dipole_matrix_config_entries(boost::program_options::options_description& desc);

	// config value getters
	const dipole_matrix_type& dipole_matrix_plus() const
	{ return dipole_matrix_plus_; }

	const dipole_matrix_type& dipole_matrix_minus() const
	{ return dipole_matrix_minus_; }

protected:
	void check(int_t system_sites);

	void post_process();

protected:
	// config values
	dipole_matrix_type dipole_matrix_plus_;
	dipole_matrix_type dipole_matrix_minus_;

};

} // namespace heom

#endif // heom_dipole_matrix_config_entries_hpp
