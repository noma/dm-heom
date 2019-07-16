// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef pathway_hpp
#define pathway_hpp

#include <iostream>
#include <map>

#include <noma/typa/typa.hpp>

#include "heom/sites_to_states.hpp"

namespace heom {

enum class pathway_dipole_sign_t
{
	minus = 0,
	plus
};

enum class pathway_dipole_side_t
{
	left = 0,
	right
};

enum class pathway_sign_t
{
	plus = 0,
	minus
};

template<size_t SIZE>
class pathway_specification
{
public:
	using pathway_array = std::array<std::pair<pathway_dipole_sign_t, pathway_dipole_side_t>, 4>;

	pathway_specification(const pathway_array& pathway, const pathway_sign_t& sign, sites_to_states_mode_t sts_mode)
	 : pathway_(pathway), sign_(sign), sites_to_states_mode_(sts_mode)
	{ }

	const pathway_array& pathway() const { return pathway_; }
	const pathway_sign_t& sign() const { return sign_; }
	const sites_to_states_mode_t& sites_to_states_mode() const { return sites_to_states_mode_; }

private:
	pathway_array pathway_;
	pathway_sign_t sign_;
	sites_to_states_mode_t sites_to_states_mode_;
};

} // namespace heom

#endif // pathway_hpp
