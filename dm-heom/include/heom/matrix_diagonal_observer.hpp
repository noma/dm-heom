// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_matrix_diagonal_observer_hpp
#define heom_matrix_diagonal_observer_hpp

#include "heom/observer.hpp"

#include "heom/common.hpp"
#include "heom/instance.hpp"

namespace heom {

class matrix_diagonal_observation_view : public observation_view
{
public:
	matrix_diagonal_observation_view(const size_t dim, const real_t time, const real_t* data)
		: dim_(dim), time_(time), data_(data)
	{
	}

	virtual ~matrix_diagonal_observation_view() = default;

	virtual void print(std::ostream& os) const
	{
		auto format = real_format;
		os << format % time_;
		for (size_t i = 0; i < dim_; ++i)
			os << default_delimiter << format % data_[(dim_ + 1) * 2 * i];
	}

private:
	const size_t dim_;
	const real_t time_;
	const real_t* data_; // not owning this
};

template<typename SOLVER_T>
class matrix_diagonal_observer : public observer
{
public:
	matrix_diagonal_observer(const hierarchy_graph& hierarchy_graph, SOLVER_T& solver)
		: dim_(solver.heom_instance().states()),
		  solver_(solver)
	{
		top_buffer_ = new complex_t[solver_.heom_instance().size_hierarchy_top_byte() / sizeof(complex_t)];
	}

	virtual ~matrix_diagonal_observer()
	{
		delete[] top_buffer_;
	}

	virtual void write_header(std::ostream& os, bool newline) const
	{
		os << "time";
		for (size_t i = 0; i < dim_; ++i)
		{
			os << default_delimiter;
			os << "site_" << i;
		}
		if (newline)
			os << '\n';
	}

	virtual std::unique_ptr<observation_view> observe(const real_t time)
	{
		solver_.hierarchy_top(top_buffer_);
		return std::unique_ptr<observation_view> {new matrix_diagonal_observation_view(dim_, time, reinterpret_cast<real_t*>(top_buffer_))}; // create a streamable view
	}

private:
	const size_t dim_;
	SOLVER_T& solver_;
	complex_t* top_buffer_;
};

} // namespace heom

#endif // heom_matrix_diagonal_observer_hpp
