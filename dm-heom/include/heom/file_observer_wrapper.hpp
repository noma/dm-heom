// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_file_observer_wrapper_hpp
#define heom_file_observer_wrapper_hpp

#include "heom/observer.hpp"

#include <fstream>
#include <memory>
#include <string>

namespace heom {

/**
 * Observer implementation that wraps another observer and writes its observations to a file.
 * The wrapped observer is owned.
 */
class file_observer_wrapper : public observer
{
public:
	file_observer_wrapper(std::unique_ptr<observer>&& obs, const std::string& filename)
		: obs_(std::move(obs)), file_(filename)
	{
		if (!file_.is_open())
			throw std::runtime_error("file_observer_wrapper(): error: could not open output file: " + filename);

		write_header(file_, true);
	}

	virtual ~file_observer_wrapper() = default;

	virtual void write_header(std::ostream& os, bool newline) const
	{
		obs_->write_header(os, newline);
	}

	virtual std::unique_ptr<observation_view> observe(const real_t time)
	{
		auto ret = obs_->observe(time); // get a view
		file_ << *ret << '\n'; // stream to file
		return ret; // return view for further use
	}

private:
	std::unique_ptr<observer> obs_;
	std::ofstream file_;
};

} // namespace heom

#endif // heom_file_observer_wrapper
