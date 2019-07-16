// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/distributed_command_line.hpp"

#include <iostream>

#include "heom/common.hpp"
#include "heom/command_line_error.hpp"

namespace heom {

namespace bpo = ::boost::program_options;

distributed_command_line::distributed_command_line(int argc, char** argv, bool parse_distributed_command_line)
	: command_line(argc, argv, false)
{
	// positional options names (avoids redundancy below)
	const std::string partition_mapping_filename_name { "partition_mapping_filename" };

	// initialise program options description
	desc_.add_options()
		(partition_mapping_filename_name.c_str(), bpo::value(&partition_mapping_filename_), "Optional partition mapping file for distributed execution (also 3rd positional option)")
		;

	// initialise program positional options description
	pos_desc_.add(partition_mapping_filename_name.c_str(), 1);

	if (parse_distributed_command_line)
		parse(argc, argv);
}

void distributed_command_line::check()
{
	command_line::check(); // call direct super-class' check()

	// TODO: check if files exist
}

} // namespace heom
