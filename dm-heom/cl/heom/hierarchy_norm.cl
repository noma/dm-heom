// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
// Copyright (c) 2016-2017 Lucas Deecke, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "common.cl"

/**
 * Computes the norm of a complex matrix by squaring its elements and adding
 * them up.
 */
__kernel void hierarchy_norm(
	const real_t threshold,                 // threshold for filter
	__global real_t* restrict norms,        // output: array of norms
	__global mask_t* restrict filters,      // output: array of filters
	__global const complex_t* restrict A    // input: array of sigma-matrices
) 
{
	const int sigma_id = get_global_id(1) * get_global_size(0) + get_global_id(0);

	// skip padded work-items
	if (sigma_id >= NUM_MATRICES)
		return;

	const int offset = sigma_id * NUM_STATES * NUM_STATES;

	// reset norm, do not reset filter
	norms[sigma_id] = 0.0;

	for (int i = 0; i < NUM_STATES; ++i) // row
	{
		for (int j = 0; j < NUM_STATES; ++j) // column
		{
			const int element_id = i * NUM_STATES + j;
			norms[sigma_id] += cnorm(A[offset + element_id]);
		}
	}

	// deactivate unprotected matrices that fall below the threshold
	if (norms[sigma_id] - threshold <= 0.0 && filters[sigma_id] != MASK_PROTECTED)
		filters[sigma_id] = MASK_FALSE;
}

