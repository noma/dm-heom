// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_constants_hpp
#define heom_constants_hpp

#include "heom/common.hpp"

#include <cmath>

// math constants
// NOTE: do not touch values
namespace heom {
namespace constants {

constexpr real_t pi = M_PI;

constexpr real_t h_bar = 1.05457148e-34;
constexpr real_t k_b = 1.380650239684840458953781387048428305608616777679979316e-23;
constexpr real_t invcm_to_joule = 1.986456217327470903861498336971100097137908380913187891e-23;
constexpr real_t debye_to_coulomb_meter = 3.3356409519815204e-30;

constexpr real_t tss_ground_state_energy = 100000; // used to avoid population of the ground state

// constants in this file generated for bosonic pade expansion following
// Wilkins 2012, Mathematica pade_expansion_0.nb
constexpr size_t max_pade = 10;
extern const long double eta_pade[];
extern const long double xi_pade[];

} // namespace constants
} // heom constants

#endif // heom_constants_hpp
