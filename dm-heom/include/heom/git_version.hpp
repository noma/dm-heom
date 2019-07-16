// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef git_version_hpp
#define git_version_hpp

#include <ostream>

extern const char GIT_VERSION_STR[];
extern const char GIT_LOCAL_CHANGES_STR[];

void write_git_version(std::ostream& out);

#endif // git_version_hpp