// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/config_base.hpp"

#include <fstream>

#include <noma/typa/parser_error.hpp>
#include <noma/typa/util.hpp>

#include "heom/common.hpp"
#include "heom/config_error.hpp"

namespace heom {

namespace bpo = ::boost::program_options;

config_base::config_base(const std::string& config_file_name, bool parse_config)
	: desc_("HEOM options")
{
	std::stringstream program_task_help;
	program_task_help << "Program task, i.e. which application the config is intended for, one of: ";
	parser::write_map_value_list(program_task_names, program_task_help);

	// initialise program options description
	desc_.add_options()
		// program
		("program.task", bpo::value(&program_task_)->required(), program_task_help.str().c_str())
		;

	if (parse_config)
		parse(config_file_name, true); // NOTE: the base allows unregistered, as it is intended to check the program_task
}

void config_base::parse(const std::string& config_file_name, bool allow_unregistered)
{
	DEBUG_ONLY( std::cout << "heom::config_base::parse(): parsing HEOM configuration from file: " << config_file_name << std::endl; )
	try {
		std::ifstream config_file(config_file_name);
		if (config_file.fail())
			throw config_error("could not open config file: " + config_file_name);

		bpo::variables_map vm;
		bpo::store(bpo::parse_config_file(config_file, desc_, allow_unregistered), vm);
		bpo::notify(vm);
	} catch (const parser::parser_error& e) {
		throw config_error(e.what());
	} catch (const config_error& e) {
		throw e;
	} catch (const bpo::invalid_option_value& e) {
		throw config_error(std::string("faulty value: ") + e.what());
	} catch (const std::exception& e) {
		throw config_error(e.what());
	}

	check_final(); // see comment on method declaration
	check();
	post_process();
}

void config_base::check()
{

}

void config_base::check_final()
{

}

void config_base::post_process()
{

}

void config_base::help(std::ostream& out) const
{
	out << desc_ << std::endl;
}

} // namespace heom
