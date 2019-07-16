// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
// Copyright (c) 2016-2017 Lucas Deecke, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "common.cl"

/**
 * Compute dot product between two hierary buffers, interpreted as vectors.
 */
__kernel void hierarchy_arithmetic_dot_product(
	__global const complex_t* restrict A,    // input: hierarchy one
	__global const complex_t* restrict B,    // input: hierarchy two
	__global complex_t* restrict result      // output: dot(A, B)
) 
{
	int sigma_id = get_global_id(1) * get_global_size(0) + get_global_id(0);

	if (sigma_id >= NUM_MATRICES)
		return;

	int offset = sigma_id * NUM_STATES * NUM_STATES;

	for (int i = 0; i < NUM_STATES * NUM_STATES; ++i)
	{
		result[offset + i] = cmult(A[offset + i], conj(B[offset + i]));
	}
}
