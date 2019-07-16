// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_matrix_trace_observer_hpp
#define heom_matrix_trace_observer_hpp

#include "heom/observer.hpp"

#include <boost/numeric/ublas/matrix.hpp>

#include "heom/common.hpp"

namespace heom {

namespace blas = boost::numeric::ublas;
using matrix_t = blas::matrix<complex_t>;

class matrix_trace_observation_view : public observation_view
{
public:
	matrix_trace_observation_view() {}

	matrix_trace_observation_view(const real_t time, const complex_t& trace)
		: time_(time), trace_(trace)
	{}

	virtual ~matrix_trace_observation_view() = default;

	virtual void print(std::ostream& os) const
	{
		auto format = real_format;
		os << format % time_;
		os << default_delimiter;
		os << format % trace_.real();
		os << default_delimiter;
		os << format % trace_.imag();
	}

	// non observation_view interface
	const real_t& time() const { return time_; }
	const complex_t& trace() const { return trace_; }

	// check for negative time to see whether the value was already set or not
	bool unset() const { return time_ < 0.0; }

	// for accumulating
	matrix_trace_observation_view& operator+=(const matrix_trace_observation_view& other)
	{
		if (unset()) {
			time_ = other.time_;
		} else {
			bool time_equal = time_ == other.time_;
			if (!time_equal)
				std::cerr << "time mismatch: " << time_ << " vs. " << other.time_ << std::endl;
			assert(time_equal); // NOTE: intended floating point equality check
		}
		trace_ += other.trace_;
		return *this;
	}

	// for averaging after accumulating
	matrix_trace_observation_view& avg(int_t count) {
		trace_ = complex_t { trace_.real() / count, trace_.imag() / count };
		return *this;
	}

private:
	real_t time_ = -1.0; // negative means unset
	complex_t trace_ = 0.0;
};

template<typename SOLVER_T>
class matrix_trace_observer : public observer
{
public:
	/**
	 * Ctor without specifying a product matrix to multiply before computing
	 * the trace. Uses identity matrix.
	 */
	matrix_trace_observer(const hierarchy_graph& hierarchy_graph, SOLVER_T& solver)
		: dim_(solver.heom_instance().states()),
		  id_(true),
		  solver_(solver)
	{
		top_buffer_ = new complex_t[solver_.heom_instance().size_hierarchy_byte() / sizeof(complex_t)];
	}

	/**
	 * Ctor which specifies an additional product matrix to be multiplied from the left
	 * before computing the trace.
	 */
	matrix_trace_observer(const hierarchy_graph& hierarchy_graph, SOLVER_T& solver, const complex_t* prod_matrix_data)
		: dim_(solver.heom_instance().states()), prod_matrix_(new matrix_t(dim_, dim_)), solver_(solver)
	{
		for (size_t i = 0; i < dim_; ++i)
			for (size_t j = 0; j < dim_; ++j)
				(*prod_matrix_)(i,j) = prod_matrix_data[index(i,j)];

		top_buffer_ = new complex_t[solver_.heom_instance().size_hierarchy_byte() / sizeof(complex_t)];
	}

	virtual ~matrix_trace_observer()
	{
		delete[] top_buffer_;
	}

	virtual void write_header(std::ostream& os, bool newline) const
	{
		os << "time" << default_delimiter << "trace_real" << default_delimiter << "trace_imag";
		if (newline)
			os << '\n';
	}

	virtual std::unique_ptr<observation_view> observe(const real_t time)
	{
		solver_.hierarchy_top(top_buffer_);
		return std::unique_ptr<observation_view> { new matrix_trace_observation_view(time, compute_trace(reinterpret_cast<real_t*>(top_buffer_))) }; // create a streamable view
	}

	/**
	 * Return actual, not the interface type.
	 */
	matrix_trace_observation_view observe_trace(const real_t time)
	{
		solver_.hierarchy_top(top_buffer_);
		return matrix_trace_observation_view(time, compute_trace(reinterpret_cast<real_t*>(top_buffer_))); // directly return an object here
	}

private:
	size_t index(size_t i, size_t j) const
	{
		return i * dim_ + j;
	}

	complex_t compute_trace(const real_t* data)
	{
		// compute trace
		const complex_t* data_complex = reinterpret_cast<const complex_t*>(data);
		matrix_t rho(dim_, dim_);

		for (size_t i = 0; i < dim_; ++i)
			for (size_t j = 0; j < dim_; ++j)
				rho(i,j) = data_complex[index(i,j)];

		matrix_t tmp = id_ ? rho : blas::prod(*prod_matrix_, rho); // either use rho, or compute product
		//	blas::matrix_vector_range<matrix_t> diag(tmp, blas::range(0, system_sites), blas::range(0, system_sites));
		//	complex_t trace = blas::sum(diag);

		complex_t trace {0.0, 0.0};
		for (size_t i = 0; i < dim_; ++i)
			trace += tmp(i,i);

		return trace;
	}

	const size_t dim_;
	bool id_ = false;
	std::unique_ptr<matrix_t> prod_matrix_;
	SOLVER_T& solver_;
	complex_t* top_buffer_;
};

} // namespace heom

#endif // heom_matrix_trace_observer_hpp
