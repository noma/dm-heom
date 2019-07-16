// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
// Copyright (c) 2015 Tobias Kramer, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_integer_partition_hpp
#define heom_integer_partition_hpp

#include "heom/common.hpp"

namespace heom {
namespace ip {

int_t size_partiton_restrict(int_t n, int_t k);

int_t size_my_integer_partitions(int_t n, int_t s);

int_t my_integer_partitions(int_t n, int_t s, int_t* res);

} // namespace ip
} // namespace heom

#endif // heom_integer_partition_hpp
