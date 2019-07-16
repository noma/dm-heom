// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_two_dimensional_spectra_config_hpp
#define heom_two_dimensional_spectra_config_hpp

#include <boost/program_options.hpp>

#include "heom/common.hpp"
#include "heom/dipole_config_entries.hpp"
#include "heom/laser_config_entries.hpp"
#include "heom/population_dynamics_config.hpp"
#include "heom/spectra_pathway.hpp"

namespace heom {

/**
 * This class contains all options needed by the fluorescence application
 */
class two_dimensional_spectra_config : public config, public dipole_config_entries
{
public:
	/** Construct from file.
	 * config_file_name .. path to config file
	 * parse .. set to false to avoid automatic parsing (i.e. when called by a subclass)
	 */
	two_dimensional_spectra_config(const std::string& config_file_name, bool parse_config = true);

	virtual ~two_dimensional_spectra_config() {};

	// config value getters

	const int_t& spectra_steps_t_1() { return spectra_steps_t_1_;}

	const int_t& spectra_steps_t_3() { return spectra_steps_t_3_;}

	const int_t& spectra_steps_t_delay() { return spectra_steps_t_delay_;}

	using spectra_pathway_list = parser::vector_wrapper<spectra_pathway_t>;

	const spectra_pathway_list& spectra_pathways() const
	{ return spectra_pathways_; }

protected:
	virtual void check();

	virtual void check_final();

	virtual void post_process();

private:
	// config values
	int_t spectra_steps_t_1_;
	int_t spectra_steps_t_3_;
	int_t spectra_steps_t_delay_;

	spectra_pathway_list spectra_pathways_;
};

} // namespace heom

#endif // heom_two_dimensional_spectra_config_hpp
