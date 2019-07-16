// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_hierarchy_norm_observer_hpp
#define heom_hierarchy_norm_observer_hpp

#include "heom/observer.hpp"

#include <noma/typa/std_vector.hpp>

#include "heom/common.hpp"
#include "heom/hierarchy_graph.hpp"
#include "heom/instance.hpp"

namespace heom {

class hierarchy_norm_observation_view : public observation_view
{
public:
	hierarchy_norm_observation_view(const size_t dim, const real_t time, const real_t* data)
		: dim_(dim), time_(time), data_(data)
	{
	}

	virtual ~hierarchy_norm_observation_view() = default;

	virtual void print(std::ostream& os) const
	{
		auto format = real_format;
		os << format % time_;
		for (size_t i = 0; i < dim_; ++i)
			os << default_delimiter << format % data_[i];
	}

private:
	const size_t dim_; // number of hierarchy nodes, i.e. number of matrices
	const real_t time_;
	const real_t* data_; // not owning this
};

template<typename SOLVER_T>
class hierarchy_norm_observer : public observer
{
public:
	hierarchy_norm_observer(const hierarchy_graph& hierarchy_graph, SOLVER_T& solver)
		: dim_(solver.heom_instance().matrices()),
		  hierarchy_graph_(hierarchy_graph),
		  solver_(solver)
	{
		norm_buffer_ = new real_t[dim_];
	}

	virtual ~hierarchy_norm_observer()
	{
		delete [] norm_buffer_;
	}

	virtual void write_header(std::ostream& os, bool newline) const
	{
		using parser::operator<<;
		os << "time";
		for (auto& tuple : hierarchy_graph_.tuples())
			os << default_delimiter << tuple;
		if (newline)
			os << '\n';
	}

	virtual std::unique_ptr<observation_view> observe(const real_t time)
	{
		solver_.update_hierarchy_norm();
		solver_.hierarchy_norm(norm_buffer_);
		return std::unique_ptr<observation_view> { new hierarchy_norm_observation_view(dim_, time, norm_buffer_) }; // create a streamable view
	}

private:
	const size_t dim_; // number of hierarchy nodes, i.e. number of matrices
	const hierarchy_graph& hierarchy_graph_;
	SOLVER_T& solver_;
	real_t* norm_buffer_;
};

} // namespace heom

#endif // heom_hierarchy_norm_observer_hpp
