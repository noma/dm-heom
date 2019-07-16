// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include <exception>

#ifndef heom_handle_main_exception_hpp
#define heom_handle_main_exception_hpp

namespace heom {

/**
 * Intended to handle any uncaught exception at the end of main.
 *
 * Usage:
 * std::exception_ptr eptr;
 * try {
 *  // main code
 * } catch (...) {
 *  eptr = std::current_exception();
 * }
 * return handle_main_exception(eptr);
 */
int handle_main_exception(std::exception_ptr eptr);

} // namespace heom

#endif // heom_handle_main_exception_hpp