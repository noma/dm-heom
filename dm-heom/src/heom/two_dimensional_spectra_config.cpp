// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/two_dimensional_spectra_config.hpp"

namespace heom {

namespace bpo = ::boost::program_options;

two_dimensional_spectra_config::two_dimensional_spectra_config(const std::string& config_file_name, bool parse_config)
	: config(config_file_name, false),
	  dipole_config_entries(desc_)
{
	std::stringstream spectra_pathways_help;
	spectra_pathways_help << "List of pathways to compute, format: {value, ...}. Valid values are: ";
	parser::write_map_value_list(spectra_pathway_names, spectra_pathways_help);


	// add options to boost program options description
	desc_.add_options()
		// spectra
		("spectra.steps_t_1", bpo::value(&spectra_steps_t_1_)->required(), "Number of steps for TODO.") // TODO: description
		("spectra.steps_t_3", bpo::value(&spectra_steps_t_3_)->required(), "Number of steps for TODO.") // TODO: description
		("spectra.steps_t_delay", bpo::value(&spectra_steps_t_delay_)->required(), "Number of steps for TODO.") // TODO: description
		("spectra.pathways", bpo::value(&spectra_pathways_)->required(), spectra_pathways_help.str().c_str())
	;

	if (parse_config)
		parse(config_file_name);
}

void two_dimensional_spectra_config::check()
{
	config::check(); // call direct config super-class' check()

	// call generic checks from entry classes
	dipole_config_entries::check(system_sites());

	// spectra_steps_t_1
	if ((spectra_steps_t_1() % program_observe_steps()) != 0)
		throw config_error("spectra.steps_t_1 must be divisible by program.oberserve_steps");
	// spectra_steps_t_3
	if ((spectra_steps_t_3() % program_observe_steps()) != 0)
		throw config_error("spectra.steps_t_3 must be divisible by program.oberserve_steps");
}

void two_dimensional_spectra_config::check_final()
{
	// check program_task_
	if (program_task() != program_task_t::TWO_DIMENSIONAL_SPECTRA)
		throw config_error("wrong program_task for this application");

	// check observations specification, transient_absorption is a special case
	// as it averages multiple observations
	if (observations().get().size() > 1)
		throw config_error("two_dimensional_spectra only supports a single observation of type matrix_trace_two_dimensional_spectra");
	for (auto& obs : observations().get())
	{
		if (obs.get().first != heom::observation_type_t::matrix_trace_two_dimensional_spectra)
			throw config_error("invalid observation type found in observations, only matrix_trace_two_dimensional_spectra is allowed");
	}
}

void two_dimensional_spectra_config::post_process()
{
	config::post_process(); // call direct super-class' post_process()

	// call post_process from entry classes
	dipole_config_entries::post_process();
}

} // namespace heom
