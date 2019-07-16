// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_command_line_hpp
#define heom_command_line_hpp

#include <string>

#include <boost/program_options.hpp>

namespace heom {

/**
 * HEOM command line parser.
 */
class command_line
{
public:

	/**
	 * parse .. set to false to avoid automatic parsing (i.e. when called by a subclass)
	 */
	command_line(int argc, char** argv, bool parse_command_line = true);

	virtual ~command_line() {};

	/**
	 * Output help text.
	 */
	void print_help(std::ostream& out) const;

	// config value getters
	const bool& help() const { return help_; }
	const bool& no_progress() const { return no_progress_; }

	const std::string& ocl_config_filename() const { return ocl_config_filename_; }
	const std::string& heom_config_filename() const { return heom_config_filename_; }

protected:
	/**
	 * Helper to throw an exception with usage message.
	 */
	void throw_usage(const std::string& msg);

	/**
	 * Trigger boost program options parsing using the specified command line options from main().
	 * Intended to be called by a ctor.
	 */
	void parse(int argc, char** argv);

	/**
	 * Check parsed values for correctness/suitability/consistenty/...
	 * virtual to make sure super-class parse() calls sub-class checks
	 */
	virtual void check();

	boost::program_options::options_description desc_;
	boost::program_options::positional_options_description pos_desc_;
	boost::program_options::variables_map vm_;

private:
	// command line parameter options
	bool help_ = false;
	bool no_progress_ = false;

	// TODO: add stuff, e.g. checkpointing (interval, reload previous, etc.), walltime for the job (e.g. to abort early if too small, or to checkpoint accordingly...)

	// positional options
	std::string ocl_config_filename_;
	std::string heom_config_filename_;
};

} // namespace heom

#endif // heom_command_line_hpp
