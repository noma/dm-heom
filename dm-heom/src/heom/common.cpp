// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/common.hpp"

#include <cstring> // memcpy()
#include <cmath> // abs()
#include <iostream>
#include <unistd.h> // isatty()

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/version.hpp>

#include "heom/compiler_version.hpp"
#include "heom/git_version.hpp"

namespace heom {

void write_to_file(const std::string& filename, std::function<void(std::ostream&)> write_func)
{
	std::ofstream file{filename};
	if (!file.is_open()) {
		throw std::runtime_error("Error: could not open output file: " + filename);
	} else {
		write_func(file);
	}
}

void write_time_statistics_summary(const bmt::statistics& stats, const std::string& name, std::ostream& os, size_t max_name_width)
{
	size_t space_count = max_name_width > name.length() ? max_name_width - name.length() : 0;
	os << "time in " << name << ": " << std::string(space_count, ' ')
	   << boost::format("%11.2f") % std::chrono::duration_cast<bmt::seconds>(stats.sum()).count() << " s "
	   << '('
	   << "count: " << boost::format("%6i") % stats.count()
	   << ", "
	   << "average: "
	   << boost::format("%10.2f") % std::chrono::duration_cast<bmt::milliseconds>(stats.average()).count() << " ms"
	   << ')' << std::endl;
}

void write_runtime_summary_table(const std::string& title, const std::vector<runtime_summary_table_entry>& table_entries, std::ostream& os)
{
	// header
	const char col_delimiter = '|';
	const char row_delimiter = '-';
	const size_t name_col_length = 30;
	const size_t count_col_length = 9;
	const size_t sum_best_col_length = 14;
	const size_t sum_worst_col_length = 14;
	const size_t avg_best_col_length = 14;
	const size_t avg_worst_col_length = 14;

	const size_t line_width = name_col_length + 1
	                          + count_col_length + 1
	                          + sum_best_col_length + 1
	                          + sum_worst_col_length + 1
	                          + avg_best_col_length + 1
	                          + avg_worst_col_length;

	const size_t title_space = (line_width - title.length()) / 2;

	const std::string best_rank_str = "rank " + boost::lexical_cast<std::string>(table_entries[0].best().rank());
	const std::string worst_rank_str = "rank " + boost::lexical_cast<std::string>(table_entries[0].worst().rank());

	os << std::string(line_width, row_delimiter) << '\n'; // horizontal line
	os << std::string(title_space, ' ') << title << '\n'; // table name
	os << std::string(line_width, row_delimiter) << '\n'; // horizontal line

	os << std::string(name_col_length, ' ') << col_delimiter
	   << std::string(count_col_length, ' ') << col_delimiter
	   << boost::format(" %" + boost::lexical_cast<std::string>(sum_best_col_length - 2) + "s ") % "sum (best)"
	   << col_delimiter
	   << boost::format(" %" + boost::lexical_cast<std::string>(sum_worst_col_length - 2) + "s ") % "sum (worst)"
	   << col_delimiter
	   << boost::format(" %" + boost::lexical_cast<std::string>(avg_best_col_length - 2) + "s ") % "avg (best)"
	   << col_delimiter
	   << boost::format(" %" + boost::lexical_cast<std::string>(avg_worst_col_length - 1) + "s") % "avg (worst)"
	   << '\n'; // first header line

	os << boost::format("%-" + boost::lexical_cast<std::string>(name_col_length - 1) + "s ") % "timer" << col_delimiter
	   << boost::format(" %" + boost::lexical_cast<std::string>(count_col_length - 2) + "s ") % "count" << col_delimiter
	   << boost::format(" %" + boost::lexical_cast<std::string>(sum_best_col_length - 2) + "s ") % best_rank_str
	   << col_delimiter
	   << boost::format(" %" + boost::lexical_cast<std::string>(sum_worst_col_length - 2) + "s ") % worst_rank_str
	   << col_delimiter
	   << boost::format(" %" + boost::lexical_cast<std::string>(avg_best_col_length - 2) + "s ") % best_rank_str
	   << col_delimiter
	   << boost::format(" %" + boost::lexical_cast<std::string>(avg_worst_col_length - 1) + "s") % worst_rank_str
	   << '\n'; // second header line


	os << std::string(line_width, row_delimiter) << '\n'; // horizontal line

	for (auto& entry : table_entries) {
		os << boost::format("%-" + boost::lexical_cast<std::string>(name_col_length - 1) + "s ") % entry.name()
		   << col_delimiter
		   << boost::format(" %" + boost::lexical_cast<std::string>(count_col_length - 2) + "i ") % entry.best().count()
		   << col_delimiter
		   << boost::format(" %" + boost::lexical_cast<std::string>(sum_best_col_length - 4) + ".2f s ")
		      % std::chrono::duration_cast<bmt::seconds>(entry.best().sum()).count() << col_delimiter
		   << boost::format(" %" + boost::lexical_cast<std::string>(sum_worst_col_length - 4) + ".2f s ")
		      % std::chrono::duration_cast<bmt::seconds>(entry.worst().sum()).count() << col_delimiter
		   << boost::format(" %" + boost::lexical_cast<std::string>(avg_best_col_length - 5) + ".2f ms ")
		      % std::chrono::duration_cast<bmt::milliseconds>(entry.best().average()).count() << col_delimiter
		   << boost::format(" %" + boost::lexical_cast<std::string>(avg_worst_col_length - 4) + ".2f ms")
		      % std::chrono::duration_cast<bmt::milliseconds>(entry.worst().average()).count() << '\n'; // entry line
	}

	os << std::string(line_width, '-') << '\n';
}

void write_memory_summary_table(const std::string& title, const std::vector<memory_summary_table_entry>& table_entries, std::ostream& os)
{
	// header
	const char col_delimiter = '|';
	const char row_delimiter = '-';
	const size_t name_col_length = 50;
	const size_t count_col_length = 9;
	const size_t size_best_col_length = 19;
	const size_t size_worst_col_length = 19;

	const size_t line_width = name_col_length + 1
	                          + count_col_length + 1
	                          + size_best_col_length + 1
	                          + size_worst_col_length;

	const size_t title_space = (line_width - title.length()) / 2;

	const std::string best_rank_str = "rank " + boost::lexical_cast<std::string>(table_entries[0].best().rank());
	const std::string worst_rank_str = "rank " + boost::lexical_cast<std::string>(table_entries[0].worst().rank());

	os << std::string(line_width, row_delimiter) << '\n'; // horizontal line
	os << std::string(title_space, ' ') << title << '\n'; // table name
	os << std::string(line_width, row_delimiter) << '\n'; // horizontal line

	os << std::string(name_col_length, ' ') << col_delimiter
	   << std::string(count_col_length, ' ') << col_delimiter
	   << boost::format(" %" + boost::lexical_cast<std::string>(size_best_col_length - 2) + "s ") % "size (best)"
	   << col_delimiter
	   << boost::format(" %" + boost::lexical_cast<std::string>(size_worst_col_length - 1) + "s") % "size (worst)"
	   << '\n'; // first header line

	os << boost::format("%-" + boost::lexical_cast<std::string>(name_col_length - 1) + "s ") % "timer" << col_delimiter
	   << boost::format(" %" + boost::lexical_cast<std::string>(count_col_length - 2) + "s ") % "count" << col_delimiter
	   << boost::format(" %" + boost::lexical_cast<std::string>(size_best_col_length - 2) + "s ") % best_rank_str
	   << col_delimiter
	   << boost::format(" %" + boost::lexical_cast<std::string>(size_worst_col_length - 1) + "s ") % worst_rank_str
	   << '\n'; // second header line

	os << std::string(line_width, row_delimiter) << '\n'; // horizontal line

	for (auto& entry : table_entries) {
		os << boost::format("%-" + boost::lexical_cast<std::string>(name_col_length - 1) + "s ") % entry.name()
		   << col_delimiter
		   << boost::format(" %" + boost::lexical_cast<std::string>(count_col_length - 2) + "i ") % entry.best().count()
		   << col_delimiter
		   << boost::format(" %" + boost::lexical_cast<std::string>(size_best_col_length - 6) + ".2f MiB ")
		      % (entry.best().size() / 1024.0 / 1024.0) << col_delimiter
		   << boost::format(" %" + boost::lexical_cast<std::string>(size_worst_col_length - 5) + ".2f MiB ")
		      % (entry.worst().size() / 1024.0 / 1024.0) << '\n'; // entry line
	}

	os << std::string(line_width, '-') << '\n';
}

void write_compile_config(std::ostream& out)
{
	write_git_version(out);
	write_compiler_version(out);

	out << "NDEBUG                        "
	    #ifdef NDEBUG
	    << "defined"
	    #else
	    << "undefined"
	    #endif
	    << std::endl;

	out << "BOOST_VERSION                 "
	    << BOOST_VERSION
	    << std::endl;

	out << "NOMA_MEMORY_DEFAULT_ALIGNMENT "
	    << NOMA_MEMORY_DEFAULT_ALIGNMENT
	    << std::endl;

	out << "HEOM_SINGLE_PRECISION         "
	    #ifdef HEOM_SINGLE_PRECISION
	    << "defined"
	    #else
	    << "undefined"
	    #endif
	    << std::endl;
}

void write_complex_matrix(const complex_t* host_buffer, size_t size_array, std::ostream& stream)
{
	auto format = complex_format;
	const real_t* buffer = reinterpret_cast<const real_t*>(host_buffer);
	for (size_t i = 0; i < size_array * size_array; ++i) {
		stream << format % buffer[2 * i + 0] % buffer[2 * i + 1] << default_delimiter;
		if (((i + 1) % size_array) == 0)
			stream << std::endl;
	}
	stream << std::endl;
}

void write_progress(size_t current_step, size_t total_steps, const std::string& prefix, std::ostream& stream)
{
	// default for only one solver run
	write_progress(current_step, total_steps, 1, 1, prefix, stream);
}

void write_progress(size_t current_step, size_t total_steps, size_t current_run, size_t total_runs, const std::string& prefix, std::ostream& stream)
{
	const bool is_stdout_and_tty = (stream.rdbuf() == std::cout.rdbuf()) &&
	                               isatty(fileno(stdout)); // stream is std::cout or an aliased stream, and is std::cout is a tty

	double percentage_current_run = (current_step + 1) * 100.0 / total_steps;
	double percentage_overall = (current_run - 1) * 100.0 / total_runs + percentage_current_run / total_runs;

	if (is_stdout_and_tty)
		std::cout << '\r'; // overwrite last status line if we're on a terminal

	size_t total_runs_str_length = boost::lexical_cast<std::string>(total_runs).size(); // number of digits of total_runs as a string
	std::string total_runs_str_length_str = boost::lexical_cast<std::string>(total_runs_str_length); // the number of digits as string to construct the format string below for equally long outputs
	std::cout << prefix << "solver run " <<  boost::format("[%" + total_runs_str_length_str + "d/%" + total_runs_str_length_str + "d]") % current_run % total_runs;
	std::cout << " at ";
	std::cout << boost::format("[%6.2f%%]") % percentage_current_run;
	std::cout << ", overall: ";
	std::cout << boost::format("[%6.2f%%]") % percentage_overall;

	if (is_stdout_and_tty)
		std::cout << std::flush;
	else
		std::cout << std::endl;

	if (current_step + 1 == total_steps && is_stdout_and_tty)
		std::cout << std::endl;
}

void write_progress_and_convergence(size_t current_step, size_t total_steps, real_t convergence, const std::string& prefix, std::ostream& stream)
{
	const bool is_stdout_and_tty = (stream.rdbuf() == std::cout.rdbuf()) &&
	                               isatty(fileno(stdout)); // stream is std::cout or an aliased stream, and is stdcout is a tty

	double current_percentage = (current_step + 1) * 100.0 / total_steps;

	if (is_stdout_and_tty)
		std::cout << '\r'; // overwrite last status line if we're on a terminal

	std::cout << prefix
	          << boost::format("[%6.2f%%]") % current_percentage
	          << " "
	          << boost::format("[%5.2e]") % convergence;

	if (is_stdout_and_tty)
		std::cout << std::flush;
	else
		std::cout << std::endl;

	if (current_step + 1 == total_steps && is_stdout_and_tty)
		std::cout << std::endl;
}


void write_site_occupancy_header(size_t sites, std::ostream& stream)
{
	stream << "time";
	for (size_t i = 0; i < sites; ++i) {
		stream << default_delimiter;
		stream << "site_" << i;
	}
	stream << std::endl;
}

void write_site_occupancy(const complex_t* host_buffer, size_t sites, std::ostream& stream)
{
	auto format = real_format;
	const real_t* buffer = reinterpret_cast<const real_t*>(host_buffer);
	for (size_t i = 0; i < sites; ++i) {
		stream << format % buffer[(sites + 1) * 2 * i];
		if (i != sites - 1)
			stream << default_delimiter;
	}
	stream << std::endl;
}

void write_hierarchy_top_header(size_t sites, std::ostream& stream)
{
	stream << "time";
	for (size_t i = 0; i < sites; ++i)
		for (size_t j = 0; j < sites; ++j)
			for (size_t k = 0; k < 2; ++k) // complex
			{
				stream << default_delimiter;
				stream << "elem_" << i << "_" << j << (k % 2 == 0 ? "_real" : "_imag");
			}
	stream << std::endl;
}

void write_hierarchy_top(const complex_t* host_buffer, size_t states, std::ostream& stream)
{
	auto format = real_format;
	const real_t* buffer = reinterpret_cast<const real_t*>(host_buffer);
	const size_t buffer_size_real = states * states * 2;
	for (size_t i = 0; i < buffer_size_real; ++i) {
		stream << format % buffer[i];
		if (i != buffer_size_real - 1)
			stream << default_delimiter;
	}
	stream << std::endl;
}

void write_norm(const real_t* host_buffer, size_t size_array, std::ostream& stream)
{
	auto format = real_format;
	for (size_t i = 0; i < size_array; ++i) {
		stream << format % host_buffer[i];
		if (i != size_array - 1) { stream << default_delimiter; }
	}
	stream << std::endl;
}

} // namespace heom