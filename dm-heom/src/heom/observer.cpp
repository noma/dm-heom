// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/observer.hpp"

namespace heom {

std::ostream& operator<<(std::ostream& os, const observation_view& view)
{
	view.print(os);
	return os;
}

} // namespace heom
