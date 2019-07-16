// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include <exception>
#include <iostream>

#include "heom/circular_dichroism_config.hpp"
#include "heom/common.hpp"
#include "heom/handle_main_exception.hpp"
#include "heom/linear_absorption_config.hpp"
#include "heom/ocl_config.hpp"
#include "heom/population_dynamics_config.hpp"
#include "heom/static_fluorescence_config.hpp"
#include "heom/transient_absorption_config.hpp"
#include "heom/thermal_state_search_config.hpp"
#include "heom/two_dimensional_spectra_config.hpp"

template<typename CONFIG>
void print_config_help(const std::string& name)
{
	std::cout << "Entries supported by the HEOM configuration file for " << name << ":" << std::endl;
	CONFIG cfg {"", false};
	cfg.help(std::cout);
}

int main(int argc, char* argv[])
{
	std::exception_ptr eptr;
	try {
		// print syntax explanation
		std::cout << "Boost.Program_options generated help texts below are in command line syntax.\n"
		          << "The config file format is as follows:\n"
		          << "\t- one entry per line:\n"
		          << "\t\tname=value\n"
		          << "\t- comment character is '#'\n"
		          << "\t- entries with prefixes can be grouped in sections:\n"
		          << "\t\tprefix.name_a=value\n"
		          << "\t\tprefix.name_b=value\n"
		          << "\t\t# is equal to:\n"
		          << "\t\t[prefix]\n"
		          << "\t\tname_a=value\n"
		          << "\t\tname_b=value\n"
		          << std::endl;

		// print help texts
		print_config_help<heom::ocl_config>("OpenCL");
		print_config_help<heom::population_dynamics_config>(heom::program_task_names.at(heom::program_task_t::POPULATION_DYNAMICS));
		print_config_help<heom::thermal_state_search_config>(heom::program_task_names.at(heom::program_task_t::THERMAL_STATE_SEARCH));
		print_config_help<heom::linear_absorption_config>(heom::program_task_names.at(heom::program_task_t::LINEAR_ABSORPTION));
		print_config_help<heom::transient_absorption_config>(heom::program_task_names.at(heom::program_task_t::TRANSIENT_ABSORPTION));
		print_config_help<heom::circular_dichroism_config>(heom::program_task_names.at(heom::program_task_t::CIRCULAR_DICHROISM));
		print_config_help<heom::static_fluorescence_config>(heom::program_task_names.at(heom::program_task_t::STATIC_FLUORESCENCE));
		print_config_help<heom::two_dimensional_spectra_config>(heom::program_task_names.at(heom::program_task_t::TWO_DIMENSIONAL_SPECTRA));

	} catch (...) {
		eptr = std::current_exception();
	}

	return heom::handle_main_exception(eptr);
}
