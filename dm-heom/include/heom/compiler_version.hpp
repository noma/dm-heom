// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef compiler_version_hpp
#define compiler_version_hpp

#include <ostream>

extern const char COMPILER_ID_STR[];
extern const char COMPILER_VERSION_STR[];

void write_compiler_version(std::ostream& out);

#endif // compiler_version_hpp