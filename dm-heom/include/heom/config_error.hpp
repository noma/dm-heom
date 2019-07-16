// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_config_error_hpp
#define heom_config_error_hpp

#include <stdexcept>
#include <string>

namespace heom {

/**
 * Custom exception class for ocl configuration errors.
 */
class config_error : public std::runtime_error
{
public:
	config_error(const std::string& msg) : runtime_error("heom::config_error: " + msg) { };
};

} // namespace heom

#endif // heom_config_error_hpp
