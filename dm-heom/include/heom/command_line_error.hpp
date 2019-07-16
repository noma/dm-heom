// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_command_line_error_hpp
#define heom_command_line_error_hpp

#include <stdexcept>
#include <string>

namespace heom {

/**
 * Custom exception class for ocl configuration errors.
 */
class command_line_error : public std::runtime_error
{
public:
	//command_line_error(const std::string& msg) : runtime_error("heom::command_line_error: " + msg) { };
	command_line_error(const std::string& msg) : runtime_error(msg) { };
};

} // namespace heom

#endif // heom_command_line_error_hpp
