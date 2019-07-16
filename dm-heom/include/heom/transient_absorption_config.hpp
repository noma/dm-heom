// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_transient_absorption_config_hpp
#define heom_transient_absorption_config_hpp

#include <cstdint>
#include <fstream>

#include <boost/program_options.hpp>
#include <noma/typa/matrix.hpp>

#include "heom/common.hpp"
#include "heom/config.hpp"
#include "heom/dipole_config_entries.hpp"
#include "heom/pulse_shape.hpp"

namespace heom {

/**
 * This class contains all options needed by all HEOM-applications.
 * Special needs are addressed by sub-classes that add further options.
 */
class transient_absorption_config : public config, public dipole_config_entries
{
public:
	/** Construct from file.
	 * config_file_name .. path to config file
	 * parse .. set to false to avoid automatic parsing (i.e. when called by a subclass)
	 */
	transient_absorption_config(const std::string& config_file_name, bool parse_config = true);

	virtual ~transient_absorption_config() {};

	// config value getters
	const bool& spectra_include_esa() const
	{ return spectra_include_esa_; }

	const pulse_shape_t& probe_pulse_shape() const
	{ return probe_pulse_shape_; }

	const real_t& probe_pulse_efield() const
	{ return probe_pulse_efield_; }

	const real_t& probe_pulse_frequency() const
	{ return probe_pulse_frequency_; }

	const real_vector_t& probe_pulse_phases() const
	{ return probe_pulse_phases_; }

	const real_t& probe_pulse_time_center() const
	{ return probe_pulse_time_center_; }

	const real_t& probe_pulse_time_full_width() const
	{ return probe_pulse_time_full_width_; }

	const pulse_shape_t& pump_pulse_shape() const
	{ return pump_pulse_shape_; }

	const real_t& pump_pulse_efield() const
	{ return pump_pulse_efield_; }

	const real_t& pump_pulse_frequency() const
	{ return pump_pulse_frequency_; }

	const real_vector_t& pump_pulse_phases() const
	{ return pump_pulse_phases_; }

	const real_t& pump_pulse_time_center() const
	{ return pump_pulse_time_center_; }

	const real_t& pump_pulse_time_full_width() const
	{ return pump_pulse_time_full_width_; }

protected:
	virtual void check();

	virtual void check_final();

	virtual void post_process();

private:
	// config values
	bool spectra_include_esa_;

	pulse_shape_t probe_pulse_shape_;
	real_t probe_pulse_efield_;
	real_t probe_pulse_frequency_;
	real_vector_t probe_pulse_phases_;
	real_t probe_pulse_time_center_;
	real_t probe_pulse_time_full_width_;

	pulse_shape_t pump_pulse_shape_;
	real_t pump_pulse_efield_;
	real_t pump_pulse_frequency_;
	real_vector_t pump_pulse_phases_;
	real_t pump_pulse_time_center_;
	real_t pump_pulse_time_full_width_;
};

} // namespace heom

#endif // heom_transient_absorption_config_hpp
