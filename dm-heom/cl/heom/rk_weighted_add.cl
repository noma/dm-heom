// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "../../../thirdparty/num/cl/types.cl"

// expected defines: NUM_MATRICES, NUM_STATES

/**
 * computes the weighted sum:
 * out = y_n + h * sum(i)(coeffs(i) * k(i))
 * used in RK solver at two points:
 * k(i) = f(t + ..., HERE)
 * y(n+1) = HERE
 */
__kernel void rk_weighted_add(
	const int    n,                         // number of terms to add, i.e. coefficients to use
	const real_t h,                         // step width
	__global       real_t* restrict out,    // output buffer (accumulated on)
	__global const real_t* restrict y_n,    // input buffers 
	real_t c1,
	__global const real_t* restrict k1,     // ks as in Runge Kutta methods
	real_t c2,
	__global const real_t* restrict k2,
	real_t c3,
	__global const real_t* restrict k3,
	real_t c4,
	__global const real_t* restrict k4, 
	real_t c5,
	__global const real_t* restrict k5, 
	real_t c6,
	__global const real_t* restrict k6,
	real_t c7,
	__global const real_t* restrict k7      // TODO: add more if needed

)
{

	// kernel setup
	#define HEOM_WORK_ITEM_GRANULARITY_MATRIX
	// own uid
	//#define sigma_id (get_global_id(0))
	#define sigma_id (get_global_id(1) * get_global_size(0) + get_global_id(0))
	#define sigma_offset (sigma_id * NUM_STATES * NUM_STATES * 2)
	#define sigma_real(i, j) (sigma_offset + 2 * ((i) * NUM_STATES + (j)))
	#define sigma_imag(i, j) (sigma_offset + 2 * ((i) * NUM_STATES + (j)) + 1)
	// arbitrary uid
	#define sigma_offset_uid(u) (u * NUM_STATES * NUM_STATES * 2)
	#define sigma_uid_real(u, i, j) ((sigma_offset_uid(u)) + 2 * ((i) * NUM_STATES + (j)))
	#define sigma_uid_imag(u, i, j) ((sigma_offset_uid(u)) + 2 * ((i) * NUM_STATES + (j)) + 1)
	// hamiltonian
	#define ham_real(i, j) (2 * ((i) * NUM_STATES + (j)))
	#define ham_imag(i, j) (2 * ((i) * NUM_STATES + (j)) + 1)

	const real_t c[] = { c1, c2, c3, c4, c5, c6, c7 }; // TODO: add more if needed
	__global const real_t* restrict k[] = { k1, k2, k3, k4, k5, k6, k7 }; // TODO: add more if needed

	// process matrix elements
	int i;
	int j;

#if defined(HEOM_WORK_ITEM_GRANULARITY_MATRIX)
	for (i = 0; i < NUM_STATES; ++i) // row
	{
		for (j = 0; j < NUM_STATES; ++j) // column
		{
#elif defined(HEOM_WORK_ITEM_GRANULARITY_ELEMENT)
			i = get_global_id(1);
			j = get_global_id(2);
#else
	// fail
#endif
//			out[sigma_real(i,j)] = y_n[sigma_real(i,j)];
//			out[sigma_imag(i,j)] = y_n[sigma_imag(i,j)];
			out[sigma_real(i,j)] = 0.0;
			out[sigma_imag(i,j)] = 0.0;

			// iterate through coefficients and ks and accumulate them on out
			for (int l = 0; l < n; ++l)
			{
//				out[sigma_real(i,j)] += h * c[l] * k[l][sigma_real(i,j)];
//				out[sigma_imag(i,j)] += h * c[l] * k[l][sigma_imag(i,j)];
				if (c[l] != 0.0) 
				{
					out[sigma_real(i, j)] += c[l] * k[l][sigma_real(i, j)];
					out[sigma_imag(i, j)] += c[l] * k[l][sigma_imag(i, j)];
				}
			}
			out[sigma_real(i, j)] *= h;
			out[sigma_imag(i, j)] *= h;
			out[sigma_real(i, j)] += y_n[sigma_real(i, j)];
			out[sigma_imag(i, j)] += y_n[sigma_imag(i, j)];

#if defined(HEOM_WORK_ITEM_GRANULARITY_MATRIX)
		}
	}
#elif defined(HEOM_WORK_ITEM_GRANULARITY_ELEMENT)
	// nothing to do
#else
	// fail
#endif

}
