// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_distributed_command_line_hpp
#define heom_distributed_command_line_hpp

#include <string>

#include <boost/program_options.hpp>

#include "heom/command_line.hpp"

namespace heom {

/**
 * HEOM command line parser for distributed apps.
 */
class distributed_command_line : public command_line
{
public:

	/**
	 * parse .. set to false to avoid automatic parsing (i.e. when called by a subclass)
	 */
	distributed_command_line(int argc, char** argv, bool parse_distributed_command_line = true);

	virtual ~distributed_command_line() {};

	// config value getters
	const std::string& partition_mapping_filename() const { return partition_mapping_filename_; }

protected:
	/**
	 * Check parsed values for correctness/suitability/consistenty/...
	 * virtual to make sure super-class parse() calls sub-class checks
	 */
	virtual void check();

private:
	// positional options
	std::string partition_mapping_filename_;
};

} // namespace heom

#endif // heom_distributed_command_line_hpp
