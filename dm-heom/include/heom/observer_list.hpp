// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_observer_list_hpp
#define heom_observer_list_hpp

#include <memory>
#include <vector>

#include <noma/typa/typa.hpp>

#include "heom/observation_type.hpp"
#include "heom/observer.hpp"

namespace heom {

/**
 * List of observers that multiplexes an observation input to its elements.
 * It has an observe() method like the observer interface, but does not return an observation_view,
 * as it is undefined how to combine the observations of its members
 * The contained observers are owned.
 */
class observer_list
{
public:
	observer_list() = default;

	observer_list(size_t capacity)
	{
		observers_.reserve(capacity);
	}

	void add(std::unique_ptr<observer>&& obs)
	{
		observers_.emplace_back(std::move(obs));
	}

	void observe(const real_t time)
	{
		for (auto& obs : observers_)
			obs->observe(time);
	}

private:
	std::vector<std::unique_ptr<observer>> observers_;
};

} // namespace heom

#endif // heom_observer_list_hpp
