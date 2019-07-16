// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include <exception>
#include <iostream>

#include <boost/format.hpp>
#include <heom/static_fluorescence_config.hpp>

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

int main(int argc, char* argv[])
{
	std::string heom_config_file = "config.cfg";
	if (argc >= 2)
		heom_config_file = argv[1];

	std::exception_ptr eptr;
	try {
		// create a config_base to check the program_task
		heom::config_base config_base(heom_config_file, true);

		std::cout << "\nChecking config file '" << heom_config_file << "' for program_task '" << config_base.program_task() << "'\n" << std::endl;

		// instantiate config according to program_task, which will throw exceptions with
		// error messages if configuration files are malformed
		switch (config_base.program_task()) {
			case heom::program_task_t::TEST:
				std::cout << "Nothing to check for test configurations." << std::endl;
				break;
			case heom::program_task_t::POPULATION_DYNAMICS:
				heom::population_dynamics_config { heom_config_file };
				break;
			case heom::program_task_t::THERMAL_STATE_SEARCH:
				heom::thermal_state_search_config { heom_config_file };
				break;
			case heom::program_task_t::LINEAR_ABSORPTION:
				heom::linear_absorption_config { heom_config_file };
				break;
			case heom::program_task_t::TRANSIENT_ABSORPTION:
				heom::transient_absorption_config { heom_config_file };
				break;
			case heom::program_task_t::CIRCULAR_DICHROISM:
				heom::circular_dichroism_config { heom_config_file };
				break;
			case heom::program_task_t::STATIC_FLUORESCENCE:
				heom::static_fluorescence_config { heom_config_file };
				break;
			case heom::program_task_t::TWO_DIMENSIONAL_SPECTRA:
				heom::two_dimensional_spectra_config { heom_config_file };
				break;
			default:
				throw std::runtime_error("Error: invalid program_task encountered.");
		}

		std::cout << "No errors found in configuration." << std::endl;

	} catch (...) {
		eptr = std::current_exception();
	}

	return heom::handle_main_exception(eptr);
}
