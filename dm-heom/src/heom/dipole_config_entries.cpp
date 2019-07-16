// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/dipole_config_entries.hpp"

#include "heom/constants.hpp"
#include "heom/config_error.hpp"

namespace heom {

namespace bpo = ::boost::program_options;

dipole_config_entries::dipole_config_entries(boost::program_options::options_description& desc)
{
	// add options to boost program options description
	desc.add_options()
		("dipole.directions", bpo::value(&dipole_directions_)->required(), "List of dipole directions (Cartesian coordinates).")
		("dipole.strengths", bpo::value(&dipole_strengths_)->required(), "List of dipole strengths for each direction (in Debye).")
		("dipole.centers", bpo::value(&dipole_centers_)->required(), "List of dipole centers (Cartesian coordinates).")
		("dipole.tensor_prefactors", bpo::value(&dipole_tensor_prefactors_)->required(), "List of dipole prefactors for tensor averaging.")
		("dipole.tensor_components", bpo::value(&dipole_tensor_components_)->required(), "List of dipole.direction components for tensor averaging.")
		;
}

void dipole_config_entries::check(int_t system_sites)
{
	// TODO: states vs. sites in 2d spectra case
//	if (dipole_directions_.size() != static_cast<size_t>(system_sites))
//		throw config_error("dipole.directions must have system.sites entries.");
//
//	if (dipole_centers_.size() != static_cast<size_t>(system_sites))
//		throw config_error("dipole.centers must have system.sites entries.");

	if (dipole_directions().size() != dipole_strengths().size() || dipole_directions().size() != dipole_centers().size())
		throw config_error("dipole.directions, dipole.strengths and dipole.centers must have the same size");

	if (dipole_tensor_prefactors_.size() != dipole_tensor_components_.size())
		throw config_error("dipole.tensor_prefactors and dipole.tensor_components must have the same size");
}

void dipole_config_entries::post_process()
{
	//dipole_strengths_.scale(heom::constants::debye_to_coulomb_meter);
}

} // namespace heom
