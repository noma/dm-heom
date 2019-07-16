// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/handle_main_exception.hpp"

#include <iostream>
#include "heom/config.hpp"

namespace heom {

int handle_main_exception(std::exception_ptr eptr)
{
	int ret = 0;
	try {
		if (eptr) {
			std::rethrow_exception(eptr);
		}
	} catch (const heom::config_error& e) {
		std::cerr << e.what() << std::endl;
		ret = -1;
	} catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		ret = -1;
	} catch (...) {
		std::cerr << "Caught unhandled exception." << std::endl;
		ret = -1;
	}
	return ret;
}

} // namespace heom