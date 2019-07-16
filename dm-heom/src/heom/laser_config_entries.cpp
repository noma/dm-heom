// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/laser_config_entries.hpp"

#include "heom/constants.hpp"
#include "heom/config_error.hpp"

namespace heom {

namespace bpo = ::boost::program_options;

laser_config_entries::laser_config_entries(boost::program_options::options_description& desc)
{
	// add options to boost program options description
	desc.add_options()
		("laser.directions", bpo::value(&laser_directions_)->required(), "List of laser directions (Cartesian coordinates, normalised).")
		;
}

void laser_config_entries::check()
{
	// TODO: it's every direction with every polarization
//	if (laser_directions_.size() != laser_polarizations_.size())
//		throw config_error("laser.directions and laser.polarizations must have the same number of entries.");
}

void laser_config_entries::post_process()
{
}

} // namespace heom
