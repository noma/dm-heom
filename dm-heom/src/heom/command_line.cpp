// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/command_line.hpp"

#include <iostream>

#include "heom/command_line_error.hpp"
#include "heom/common.hpp"

namespace heom {

namespace bpo = ::boost::program_options;

command_line::command_line(int argc, char** argv, bool parse_command_line)
	: desc_("Available command line options")
{
	// positional options names (avoids redundancy below)
	const std::string ocl_config_filename_name { "ocl_config_filename" };
	const std::string heom_config_filename_name { "heom_config_filename" };

	// initialise program options description
	desc_.add_options()
		("help,h", bpo::value(&help_)->zero_tokens(), "Print help messages")
		("no-progress", bpo::value(&no_progress_)->zero_tokens(), "Do not show progress output.")
		(ocl_config_filename_name.c_str(), bpo::value(&ocl_config_filename_)->required(), "OpenCL configuration file (also 1st positional option)")
		(heom_config_filename_name.c_str(), bpo::value(&heom_config_filename_)->required(), "App-specific HEOM configuration file (also 2nd positional option)")
		;

	// initialise program positional options description
	pos_desc_.add(ocl_config_filename_name.c_str(), 1);
	pos_desc_.add(heom_config_filename_name.c_str(), 1);

	if (parse_command_line)
		parse(argc, argv);
}

void command_line::throw_usage(const std::string& msg)
{
	std::stringstream usage;
	usage << std::endl;
	print_help(usage);
	usage << "heom::command_line_error: " << msg;
	throw command_line_error(usage.str());
}

void command_line::parse(int argc, char** argv)
{
	DEBUG_ONLY( std::cout << "heom::command_line::parse(): parsing HEOM command line" << std::endl; )

	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc_).positional(pos_desc_).run(), vm_);

		if (vm_.count("help"))
			throw command_line_error("Help option was specified, just showing usage.");

		bpo::notify(vm_);

		check();
	} catch (const std::exception& e) {
		throw_usage(e.what());
	}
}

void command_line::check()
{
	// TODO: check if files exist
}

void command_line::print_help(std::ostream& out) const
{
	out << "Usage: <app_name> [options] ";
	for (size_t i = 0; i < pos_desc_.max_total_count(); ++i )
		out << '<' << pos_desc_.name_for_position(i) << '>' << ' ';
	out << std::endl << std::endl;
	out << desc_ << std::endl;
}

} // namespace heom
