// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/dipole_matrix_config_entries.hpp"

#include "heom/constants.hpp"
#include "heom/config_error.hpp"

namespace heom {

namespace bpo = ::boost::program_options;

dipole_matrix_config_entries::dipole_matrix_config_entries(boost::program_options::options_description& desc)
{
	// add options to boost program options description
	desc.add_options()
		("dipole_matrix.plus", bpo::value(&dipole_matrix_plus_)->required(), "List of dipole matrices during propagation (system_sites x system_sites).")
		("dipole_matrix.minus", bpo::value(&dipole_matrix_minus_)->required(), "List of dipole matrice for trace operation (system_sites x system_sites).")
		;
}

void dipole_matrix_config_entries::check(int_t system_sites)
{
	// dipole_matrix_plus
	for (size_t i = 0; i < dipole_matrix_plus_.size(); ++i)
		if (!(static_cast<int_t>(dipole_matrix_plus_[i].rows()) == system_sites &&
		      static_cast<int_t>(dipole_matrix_plus_[i].cols()) == system_sites))
			throw config_error("dipole_matrix.plus elements must be system.sites by system.sites matrices");

	// dipole_matrix_minus
	for (size_t i = 0; i < dipole_matrix_minus_.size(); ++i)
		if (!(static_cast<int_t>(dipole_matrix_minus_[i].rows()) == system_sites &&
		      static_cast<int_t>(dipole_matrix_minus_[i].cols()) == system_sites))
			throw config_error("dipole_matrix.minus elements must be system.sites by system.sites matrices");

	if (dipole_matrix_plus_.size() != dipole_matrix_minus_.size())
		throw config_error("dipole_matrix.plus and dipole_matrix.minus must have the same number of elements");
}

void dipole_matrix_config_entries::post_process()
{
	// scale to SI units
	for (size_t i = 0; i < dipole_matrix_plus_.size(); ++i) {
		dipole_matrix_plus_[i].scale(constants::invcm_to_joule);
		dipole_matrix_minus_[i].scale(constants::invcm_to_joule);
	}

	// TODO: maybe implement begin/end etc. for parser::vector for this to work (plus and minus)
	//for (auto& m : dipole_matrix_minus_)
	//	m.scale(constants::invcm_to_joule);
}

} // namespace heom
