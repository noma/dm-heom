// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
// Copyright (c) 2016-2017 Lucas Deecke, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "common.cl"

__kernel void error_norm(
	__global real_t* restrict error_norms,    // output: array of error norms
	__global const complex_t* restrict A,     // input: array of sigma-matrices
	__global const complex_t* restrict B      // input: array of sigma-matrices
) 
{
	const int sigma_id = get_global_id(1) * get_global_size(0) + get_global_id(0);

	// skip padded work-items
	if (sigma_id >= NUM_MATRICES)
		return;
	
	const int offset = sigma_id * NUM_STATES * NUM_STATES;

	// reset error_norm
	real_t norm = 0.0;
	error_norms[sigma_id] = 0.0;

	for (int i = 0; i < NUM_STATES; ++i) // row
	{
		for (int j = 0; j < NUM_STATES; ++j) // column
		{
			const int element_id = i * NUM_STATES + j;
			complex_t d = A[offset + element_id] - B[offset + element_id];
			norm += real(cmult(conj(d), d));
		}
	}

	error_norms[sigma_id] = sqrt(norm);
}

