// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_observer_hpp
#define heom_observer_hpp

#include <ostream>
#include <memory>

#include "heom/common.hpp"

namespace heom {

/**
 * Observation view interface.
 */
class observation_view
{
public:
	virtual ~observation_view() = default;

	virtual void print(std::ostream& os) const = 0;
};


/**
 * Observer interface.
 */
class observer
{
public:
	virtual ~observer() = default;


	virtual void write_header(std::ostream& os, bool newline) const = 0;

	/**
	 * Observe by returning a streamable observation_view object (managed by unique_ptr).
	 */
	virtual std::unique_ptr<observation_view> observe(const real_t time) = 0;
};


/**
 * Generic observation_view output stream operator, calling the view's print() member funtion.
 */
std::ostream& operator<<(std::ostream& os, const observation_view& view);

} // namespace heom

#endif // heom_observer_hpp
