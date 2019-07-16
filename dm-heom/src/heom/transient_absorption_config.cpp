// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/transient_absorption_config.hpp"

#include <fstream>
#include <iostream>

namespace heom {

namespace bpo = ::boost::program_options;

transient_absorption_config::transient_absorption_config(const std::string& config_file_name, bool parse_config)
	: config(config_file_name, false),
	  dipole_config_entries(desc_)
{
	// add options to boost program options description
	desc_.add_options()
		// interaction with a time-dependent laserfield consisting of a probe pulse
		("spectra.include_esa", bpo::value(&spectra_include_esa_)->required(), "Enable Excited State Absorption pathway.")
		("probe_pulse.shape", bpo::value(&probe_pulse_shape_)->required(), "Shape function of probe pulse [Gaussian].") // string {"Gaussian"}
		("probe_pulse.efield", bpo::value(&probe_pulse_efield_)->required(), "Peak electric field of probe pulse [V/m].")
		("probe_pulse.frequency", bpo::value(&probe_pulse_frequency_)->required(), "Frequency of probe pulse [invcm].")
		("probe_pulse.phases", bpo::value(&probe_pulse_phases_)->required(), "Phases of probe pulse [radians].")
		("probe_pulse.time_center", bpo::value(&probe_pulse_time_center_)->required(), "Center time of probe pulse [s].")
		("probe_pulse.time_full_width", bpo::value(&probe_pulse_time_full_width_)->required(), "Full width of probe pulse [s].")
		// interaction with a time-dependent laserfield consisting of a pump pulse
		("pump_pulse.shape", bpo::value(&pump_pulse_shape_)->required(), "Shape function of pump pulse [Gaussian].")
		("pump_pulse.efield", bpo::value(&pump_pulse_efield_)->required(), "Peak electric field of pump pulse [V/m].")
		("pump_pulse.frequency", bpo::value(&pump_pulse_frequency_)->required(), "Frequency of pump pulse [invcm].")
		("pump_pulse.phases", bpo::value(&pump_pulse_phases_)->required(), "Phases of pump pulse [radians].")
		("pump_pulse.time_center", bpo::value(&pump_pulse_time_center_)->required(), "Center time of pump pulse [s].")
		("pump_pulse.time_full_width", bpo::value(&pump_pulse_time_full_width_)->required(), "Full width of pump pulse [s].");;
		;

	if (parse_config)
		parse(config_file_name);
}

void transient_absorption_config::check()
{
	config::check(); // call direct super-class' check()

	// call generic checks from entry classes
	dipole_config_entries::check(system_sites());

	// probe_pulse_shape
	// TODO: adapt if more shapes are needed
	if (probe_pulse_shape_ != pulse_shape_t::GAUSSIAN)
		throw config_error("probe_pulse.shape must be '" + pulse_shape_names.at(pulse_shape_t::GAUSSIAN) + "'");

	// probe_pulse_efield_

	// probe_pulse_frequency_

	// probe_pulse_phases_

	// probe_pulse_time_center_

	// probe_pulse_time_full_width_


	// TODO: adapt if more shapes are needed
	if (pump_pulse_shape_ != pulse_shape_t::GAUSSIAN)
		throw config_error("pump_pulse.shape must be '" + pulse_shape_names.at(pulse_shape_t::GAUSSIAN) + "'");

	// pump_pulse_efield_

	// pump_pulse_frequency_

	// pump_pulse_phases_
	if (!(pump_pulse_phases_.size() == probe_pulse_phases().size()))
		throw config_error("pump_pulse.phases must have the same size as probe_pulse.phases");

	// pump_pulse_time_center_

	// pump_pulse_time_full_width_
}

void transient_absorption_config::check_final()
{
	// check program_task_
	if (program_task() != program_task_t::TRANSIENT_ABSORPTION)
		throw config_error("wrong program_task for this application");

	// check observations specification, transient_absorption is a special case
	// as it averages multiple observations
//	if (observations().get().size() > 1)
//		throw config_error("transient_absorption only supports a single observation of type matrix_trace_transient_absorption");
//	for (auto& obs : observations().get())
//	{
//		if (obs.get().first != heom::observation_type_t::matrix_trace_transient_absorption)
//			throw config_error("invalid observation type found in observations, only matrix_trace_transient_absorption is allowed");
//	}
}

void transient_absorption_config::post_process()
{
	config::post_process(); // call direct super-class' post_process()

	// call post_process from entry classes
	dipole_config_entries::post_process();
}

} // namespace heom