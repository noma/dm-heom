// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_hierarchy_complete_observer_hpp
#define heom_hierarchy_complete_observer_hpp

#include "heom/observer.hpp"

#include "heom/common.hpp"
#include "heom/instance.hpp"

namespace heom {

class hierarchy_complete_observation_view : public observation_view
{
public:
	hierarchy_complete_observation_view(const size_t hierarchy_dim, const size_t matrix_dim, const real_t time, const real_t* data)
		: size_(hierarchy_dim * matrix_dim * matrix_dim * 2), time_(time), data_(data)
	{}

	virtual ~hierarchy_complete_observation_view() = default;

	virtual void print(std::ostream& os) const
	{
		auto format = real_format;
		os << format % time_;
		for (size_t i = 0; i < size_; ++i)
			os << default_delimiter << format % data_[i];
	}

private:
	const size_t size_;
	const real_t time_;
	const real_t* data_; // not owning this
};

template<typename SOLVER_T>
class hierarchy_complete_observer : public observer
{
public:
	hierarchy_complete_observer(const hierarchy_graph& hierarchy_graph, SOLVER_T& solver)
		: hierarchy_dim_(solver.heom_instance().matrices()),
		  matrix_dim_(solver.heom_instance().states()),
		  solver_(solver)
	{
		buffer_ = new complex_t[solver_.heom_instance().size_hierarchy_byte() / sizeof(complex_t)];
	}

	virtual ~hierarchy_complete_observer()
	{
		delete[] buffer_;
	}

	virtual void write_header(std::ostream& os, bool newline) const
	{
		os << "time";
		for (size_t h = 0; h < hierarchy_dim_; ++h)
			for (size_t i = 0; i < matrix_dim_; ++i)
				for (size_t j = 0; j < matrix_dim_; ++j)
					for (size_t k = 0; k < 2; ++k) // complex
					{
						os << default_delimiter;
						os << "mat_" << h << "_elem_" << i << "_" << j << (k % 2 == 0 ? "_real" : "_imag");
					}
		if (newline)
			os << '\n';
	}

	virtual std::unique_ptr<observation_view> observe(const real_t time)
	{
		solver_.hierarchy(buffer_);
		return std::unique_ptr<observation_view> { new hierarchy_complete_observation_view(hierarchy_dim_, matrix_dim_, time, reinterpret_cast<real_t*>(buffer_)) }; // create a streamable view
	}

private:
	const size_t hierarchy_dim_;
	const size_t matrix_dim_;
	SOLVER_T& solver_;
	complex_t* buffer_;
};

} // namespace heom

#endif // heom_hierarchy_complete_observer_hpp
