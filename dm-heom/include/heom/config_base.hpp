// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_config_base_hpp
#define heom_config_base_hpp

#include <iostream>
#include <string>

#include <boost/program_options.hpp>

#include "heom/program_task.hpp"

namespace heom {

/**
 * This class is the top-level base of the config class inheritance tree, it
 * contains only a single program_task entry which allows to check a given
 * configuration file for the task, to instantiate the corresponding configuration
 * type dynamically at runtime.
 */
class config_base
{
public:

	/** Construct from file.
	 * config_file_name .. path to config file
	 * parse .. set to false to avoid automatic parsing (i.e. when called by a subclass)
	 */
	config_base(const std::string& config_file_name, bool parse_config = true);

	virtual ~config_base() {};

	/**
	 * Output help text.
	 */
	void help(std::ostream& out) const;

	// config value getters

	// program
	const program_task_t& program_task() const
	{ return program_task_; }

protected:
	/**
	 * Trigger boost program options parsing using the specified config file.
	 * Intended to be called by a config class' ctor.
	 *
	 * For actual config files that correspond to applications allow_unregistered should
	 * be set to false, to enforce errors for unknown values within the configuration file.
	 * For instantiating config_base, i.e. to check for the program_task, it can is set to
	 * true, to not fail for actual configuration files.
	 */
	void parse(const std::string& config_file_name, bool allow_unregistered = false);

	/**
	 * Check parsed values for correctness/suitability/consistenty/...
	 * virtual to make sure super-class parse() calls sub-class checks
	 */
	virtual void check();

	/**
	 * Check values that depend on the instance's actual type.
	 * E.g. linear_absorption_config needs to check its program_task
	 * only if it is not inside an instance of transient_absorption_config.
	 * To be implemented only by sub-class that can be instanciated.
	 */
	virtual void check_final();

	/**
	 * Do post-processing, like unit conversions, computing derived values etc.
	 * virtual to make sure super-class parse() calls sub-class checks
	 */
	virtual void post_process();

	boost::program_options::options_description desc_;

private:
	// config values

	// program
	program_task_t program_task_;
};

} // namespace heom

#endif // heom_config_base_hpp
