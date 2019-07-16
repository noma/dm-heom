// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_common_hpp
#define heom_common_hpp

#include <functional>
#include <iostream>
#include <fstream>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <noma/bmt/bmt.hpp> // benchmark timer
#include <noma/memory/memory.hpp>
#include <noma/num/types.hpp>
#include <noma/typa/vector.hpp>
#include <noma/typa/matrix.hpp>

// forward declaration
namespace noma {
namespace typa {}
namespace ocl {}
namespace memory {}
}

namespace heom {

namespace num = ::noma::num;

using num::bool_t;
using num::int_t;
using num::uint_t;

using num::real_t;
using num::complex_t;
using num::long_real_t;
using num::long_complex_t;

using num::int_format;
using num::real_format;
using num::complex_format;
using num::default_delimiter;

namespace parser = ::noma::typa;
namespace ocl = ::noma::ocl;
namespace memory = ::noma::memory;

using real_vector_t    = parser::vector<real_t>;
using complex_vector_t = parser::vector<complex_t>;
using real_matrix_t    = parser::matrix<real_t>;
using complex_matrix_t = parser::matrix<complex_t>;

// bring default integer types in
using std::int32_t;
using std::int64_t;

// TODO: should be in communicator, but is used by some of the code below
// identifier for processses in a distributed run, should be an integral type
using rank_t = int; // MPI compatible

// vector of statistic objects references together with a name
namespace bmt = ::noma::bmt;
using statistics_vector = std::vector<std::pair<const bmt::statistics&, const std::string>>;

class runtime_summary_record
{
public:
	runtime_summary_record() : rank_(-1), count_(0), sum_(0), average_(0)
	{}

	runtime_summary_record(rank_t rank, const bmt::statistics& stats)
		: rank_(rank), count_(stats.count()), sum_(stats.sum()), average_(stats.average())
	{}

	rank_t rank() const
	{ return rank_; }

	size_t count() const
	{ return count_; }

	const bmt::duration& sum() const
	{ return sum_; }

	const bmt::duration& average() const
	{ return average_; }

private:
	rank_t rank_;
	size_t count_;
	bmt::duration sum_;
	bmt::duration average_;
};

class memory_summary_record
{
public:
	memory_summary_record() : rank_(-1), count_(0), size_(0)
	{}

	memory_summary_record(rank_t rank, size_t count, size_t size)
		: rank_(rank), count_(count), size_(size)
	{}

	rank_t rank() const
	{ return rank_; }

	size_t count() const
	{ return count_; }

	size_t size() const
	{ return size_; }

private:
	rank_t rank_;
	size_t count_;
	size_t size_;
};

template<typename T>
class summary_table_entry
{
public:
	summary_table_entry(const std::string& name, const T& best, const T& worst)
		: name_(name), best_(best), worst_(worst)
	{}

	const std::string& name() const
	{ return name_; }

	const T& best() const
	{ return best_; }

	const T& worst() const
	{ return worst_; }

private:
	const std::string name_;
	T best_;
	T worst_;
};

using runtime_summary_table_entry = summary_table_entry<runtime_summary_record>;
using memory_summary_table_entry = summary_table_entry<memory_summary_record>;

void write_runtime_summary_table(const std::string& title, const std::vector<runtime_summary_table_entry>& table_entries, std::ostream& os);

void write_memory_summary_table(const std::string& title, const std::vector<memory_summary_table_entry>& table_entries, std::ostream& os);

// easy way to write all the functions that take a stream to a file, use std::bind if write_func has more arguments
// NOTE: we could do this way easier by some std::ofstream create_ofstream(filename), and pass it this expression as argument
//       but move semantics for streams require a GCC 5 STL
void write_to_file(const std::string& filename, std::function<void(std::ostream&)> write_func);

// generic summary output for a statistics object
void write_time_statistics_summary(const bmt::statistics& stats, const std::string& name, std::ostream& os, size_t max_name_width = 25);

// print compile time flags
void write_compile_config(std::ostream& out);

// print complex square matrix
void write_complex_matrix(const complex_t* host_buffer, size_t size_array, std::ostream& stream);

void write_progress(size_t current_step, size_t total_steps, const std::string& prefix, std::ostream& stream);
void write_progress(size_t current_step, size_t total_steps, size_t current_run, size_t total_runs, const std::string& prefix, std::ostream& stream);

void write_progress_and_convergence(size_t current_step, size_t total_steps, real_t convergence, const std::string& prefix, std::ostream& stream);

void write_site_occupancy_header(size_t size_array, std::ostream& stream);

void write_site_occupancy(const complex_t* host_buffer, size_t size_array, std::ostream& stream);

void write_hierarchy_top_header(size_t states, std::ostream& stream);

void write_hierarchy_top(const complex_t* host_buffer, size_t states, std::ostream& stream);

void write_norm(const real_t* host_buffer, size_t size_array, std::ostream& stream);

} // namespace heom

#endif // common_hpp
