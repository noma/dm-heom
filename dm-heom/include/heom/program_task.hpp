// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_program_task_hpp
#define heom_program_task_hpp

#include <iostream>
#include <map>

namespace heom {

enum class program_task_t
{
	TEST = 0,
	POPULATION_DYNAMICS,
	RATE_KERNEL,
	THERMAL_STATE_SEARCH,
	LINEAR_ABSORPTION,
	TRANSIENT_ABSORPTION,
	CIRCULAR_DICHROISM,
	STATIC_FLUORESCENCE,
	TWO_DIMENSIONAL_SPECTRA
};

const std::map<program_task_t, std::string> program_task_names { 
	{ program_task_t::TEST, "test" },
	{ program_task_t::POPULATION_DYNAMICS, "population_dynamics" },
	{ program_task_t::RATE_KERNEL, "rate_kernel" },
	{ program_task_t::THERMAL_STATE_SEARCH, "thermal_state_search" },
	{ program_task_t::LINEAR_ABSORPTION, "linear_absorption" },
	{ program_task_t::TRANSIENT_ABSORPTION, "transient_absorption" },
	{ program_task_t::CIRCULAR_DICHROISM, "circular_dichroism" },
	{ program_task_t::STATIC_FLUORESCENCE, "static_fluorescence" },
	{ program_task_t::TWO_DIMENSIONAL_SPECTRA, "two_dimensional_spectra" }
};

std::ostream& operator<<(std::ostream& out, const program_task_t& t);
std::istream& operator>>(std::istream& in, program_task_t& t);

} // namespace heom

#endif // heom_program_task_hpp
