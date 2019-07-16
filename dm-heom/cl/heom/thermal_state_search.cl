// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "common.cl"

__kernel void thermal_state_search(
	__global complex_t* restrict out,               // output: result buffer, hierarchy state
	__global const complex_t* restrict in,          // input: initial value, hierarchy state
	const real_t epsilon_real,
	const real_t epsilon_imag,
	__global const complex_t* restrict c, // TODO: remove
	__global const complex_t* restrict c_transposed, // TODO: remove
	__global const complex_t* restrict eigenvalues,
	__global const complex_t* restrict pf_same_sum
)
{
	int sigma_id = get_global_id(1) * get_global_size(0) + get_global_id(0);

	// skip padded work-items
	if (sigma_id >= NUM_MATRICES)
		return;

	// find the ADO on which to operate
	int offset = sigma_id * NUM_STATES * NUM_STATES;

	// TODO: decide if we support inplace operation for this kernel, not needed by use-case but could reduce memory bandwidth requirements
	// "out" can point to the same memory as "in", so a temporary buffer is needed
	//complex_t out_temp[NUM_STATES * NUM_STATES];

	complex_t complex_i = (complex_t)(0.0, 1.0);
	complex_t epsilon = (complex_t)(epsilon_real, epsilon_imag);

	// matrix elements
	for (int i = 0; i < NUM_STATES; ++i) {
		for (int j = 0; j < NUM_STATES; ++j) {
			complex_t value = (complex_t)(0.0, 0.0);
			int elem_id = offset + i * NUM_STATES + j;

			// formula (21), lower part
			out[elem_id] = cdiv(in[elem_id],
			                   cadd(
			                       cmult(complex_i / CONST_SI_HBAR, cadd(eigenvalues[i], -eigenvalues[j])), // (i / h_bar) * (E_i - E_j)
			                       cadd(epsilon, -pf_same_sum[sigma_id]) // epsilon + PFSAME_u
			                   )
			               );
		}
	}
}
