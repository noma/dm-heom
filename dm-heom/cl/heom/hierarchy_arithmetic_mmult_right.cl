// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
// Copyright (c) 2016-2017 Lucas Deecke, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "common.cl"

/**
 * right multiplication of each hierarchy auxiliary density matrix with input M
 * NOTE: no restrict keyword because arguments might point to the same buffers
 */
__kernel void hierarchy_arithmetic_mmult_right(
	__global const complex_t* A,    // input: hierarchy
	__global const complex_t* M,    // input: matrix with which to multiply
	__global complex_t* result      // output: hierarchy-sized, result[sigma_id] = A[sigma_id] * M
)
{
	int sigma_id = get_global_id(1) * get_global_size(0) + get_global_id(0);

	// skip padded work-items
	if (sigma_id >= NUM_MATRICES)
		return;
	
	// find the ADO on which to operate
	int offset = sigma_id * NUM_STATES * NUM_STATES;

	// result can point to the same memory as A, so a temporary buffer is needed
	complex_t result_temp[NUM_STATES * NUM_STATES];

	for (int i = 0; i < NUM_STATES; ++i)
	{
		for (int j = 0; j < NUM_STATES; ++j)
		{
			complex_t value = (complex_t)(0.0, 0.0);
			for (int k = 0; k < NUM_STATES; ++k) // result[i,j] = sum(k=1,NUM_STATES) A[i,k] * M[k,j]
			{
				value = cadd(value, cmult(A[offset + i * NUM_STATES + k], M[k * NUM_STATES + j]));
			}
			result_temp[i * NUM_STATES + j] = value;
		}
	}

	for (int i = 0; i < NUM_STATES * NUM_STATES; ++i)
	{
		result[offset + i] = result_temp[i];
	}
}
