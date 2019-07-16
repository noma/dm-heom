// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "common.cl"

__constant const int BATHS_COUPLING[] = BATHS_COUPLING_LIST;

/**
 * This kernel implements the HEOM ODE, performs an evaluation of the time derivative.
 * It has two optional features, that can also be combined: NOMA_NUM_ODE_ACCUMULATE, and NOMA_NUM_ODE_SUBDIAGONAL
 * If neither is used, step_size, out_coeff, acc, acc_coeff, acc_init are unused.
 *
 * NOMA_NUM_ODE_ACCUMULATE integrates the weighted accumulation of the k's of a Runge-Kutta scheme into their computation, which increases memory re-use and eliminates the final weighted add..
 * NOMA_NUM_ODE_ACCUMULATE does not need out_coeff
 *
 * NOMA_NUM_ODE_SUBDIAGONAL is for Runge-Kutta schemes where only the first subdiagonal is populated (e.g. classical RK4). This allows to already compute the input for the next evaluation, eliminating a weighted add.
 * NOMA_NUM_ODE_SUBDIAGONAL does not need acc, acc_coeff, acc_init
 *
 */
__kernel void heom_ode(
	const real_t step_size,                         // time step, h
	__global real_t* restrict sigma_out,            // output: result buffer, hierarchy state, evaluated derivative
	const real_t out_coeff,                         // result coefficient for writing to out (a's in butcher tableau)
	__global real_t* restrict sigma_acc,            // accumulation buffer, set to NULL to disable accumulation
	const real_t acc_coeff,                         // result coefficient for writing to acc (b's in butcher tableau)
	const bool_t acc_init,                          // if true: initialise acc to y_n, used in the first evaluation of a step
	__global const real_t* restrict sigma_in,       // input: initial value, hierarchy state, y_n
	__global const real_t* restrict hamiltonian,    // input: hamiltonian, complex matrix
	__global const int* restrict ado_tuples,        // input: ADO tuples
	__global const int* restrict ado_plus,          // input: ADO plus indices
	__global const int* restrict ado_minus,         // input: ADO minus indices
	__global const real_t* restrict pf_same_sum,    // input: precomputed parts of the equation ...
	__global const real_t* restrict pf_same,        // input: ...
	__global const real_t* restrict pf_circ,        // input: ...
	__global const real_t* restrict cbk,            // input: ...
	__global const real_t* restrict cbk_sqrt,       // input: ...
	const int first_non_plus_id,                    // input: first sigma_id that has no plus edges in the hierarchy
	__global mask_t* restrict mask,                 // input: mask to enable/disable hierarchy nodes
	__global real_t* restrict flow_from_above,      // output: flow from above
	__global real_t* restrict flow_from_below       // output: flow from below
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

	// iterate matrices
#if defined(HEOM_WORK_ITEM_GRANULARITY_MATRIX)
	for (int i = 0; i < NUM_STATES; ++i)
	{
		for (int j = 0; j < NUM_STATES; ++j)
		{
#elif defined(HEOM_WORK_ITEM_GRANULARITY_ELEMENT)
	int i = get_local_id(1);
	int j = get_local_id(2);
#else
// TODO: fail
#endif

			// -------------------- COMMUTATOR TERM --------------------
			real_t result_commutator_real = 0.0;
			real_t result_commutator_imag = 0.0;

#ifndef HEOM_ODE_DISABLE_COMMUTATOR
			// temporaries for numerically stable accumulation
			real_t sh_prod_real;
			real_t sh_prod_imag;
			real_t hs_prod_real;
			real_t hs_prod_imag;

	#ifndef HEOM_TSS_EQUATION
			for (int k = 0; k < NUM_STATES; ++k)
			{
				// NOTE: accumulation on the same variable or directly inside
				//       the result buffer (e.g result_commutator_*) are
				//       numerically less stable and lead to non-zero imaginary
				//       parts on the diagonal

				// sigma_in[i,k] * hamiltonian[k,j]
				sh_prod_real  =  sigma_in[sigma_real(i, k)] * hamiltonian[ham_real(k, j)];
				sh_prod_real -=  sigma_in[sigma_imag(i, k)] * hamiltonian[ham_imag(k, j)];
				sh_prod_imag  =  sigma_in[sigma_real(i, k)] * hamiltonian[ham_imag(k, j)];
				sh_prod_imag +=  sigma_in[sigma_imag(i, k)] * hamiltonian[ham_real(k, j)];

				// hamiltonian[k,j] - sigma_in[i,k]
				hs_prod_real  = hamiltonian[ham_real(i, k)] *  sigma_in[sigma_real(k, j)];
				hs_prod_real -= hamiltonian[ham_imag(i, k)] *  sigma_in[sigma_imag(k, j)];
				hs_prod_imag  = hamiltonian[ham_real(i, k)] *  sigma_in[sigma_imag(k, j)];
				hs_prod_imag += hamiltonian[ham_imag(i, k)] *  sigma_in[sigma_real(k, j)];

				result_commutator_real += sh_prod_real - hs_prod_real;
				result_commutator_imag += sh_prod_imag - hs_prod_imag;
			}
			
			real_t tmp = result_commutator_real;
			result_commutator_real = -result_commutator_imag / CONST_SI_HBAR;
			result_commutator_imag = tmp / CONST_SI_HBAR;

	#else // ifdef HEOM_TSS_EQUATION
			result_commutator_real  = sigma_in[sigma_real(i, j)] * hamiltonian[ham_real(i, j)];
			result_commutator_real -= sigma_in[sigma_imag(i, j)] * hamiltonian[ham_imag(i, j)];
			result_commutator_imag  = sigma_in[sigma_real(i, j)] * hamiltonian[ham_imag(i, j)];
			result_commutator_imag += sigma_in[sigma_imag(i, j)] * hamiltonian[ham_real(i, j)];
	#endif // HEOM_TSS_EQUATION

#endif // HEOM_ODE_DISABLE_COMMUTATOR

			// -------------------- SAME TERM --------------------
			real_t result_same_real  = 0.0;
			real_t result_same_imag  = 0.0;

#ifndef HEOM_ODE_DISABLE_SAME
	#ifndef HEOM_TSS_EQUATION
			// complex product of pf_same_sum[sigma_id] and sigma_in[i,j]
			const real_t pf_same_sum_real = pf_same_sum[2 * sigma_id];
			const real_t pf_same_sum_imag = pf_same_sum[2 * sigma_id + 1];
			result_same_real = pf_same_sum_real * sigma_in[sigma_real(i, j)] - pf_same_sum_imag * sigma_in[sigma_imag(i, j)];
			result_same_imag = pf_same_sum_real * sigma_in[sigma_imag(i, j)] + pf_same_sum_imag * sigma_in[sigma_real(i, j)];
	#endif // HEOM_TSS_EQUATION
			// B sum in equation (4), (5), (6)
			if (i != j) {
				for (int b = 0; b < BATHS_MAX_PER_STATE; ++b) {
					// bath b for site i
					const int b_i = BATHS_COUPLING[i * BATHS_MAX_PER_STATE + b];

					// check for coupling
					if (b_i >= 0) {
						// complex product: -pf_same[b_i] * sigma_in[id]
						const real_t minus_pf_same_real = -pf_same[2 * b_i];
						const real_t minus_pf_same_imag = -pf_same[2 * b_i + 1];
						result_same_real += minus_pf_same_real * sigma_in[sigma_real(i, j)]
						                  - minus_pf_same_imag * sigma_in[sigma_imag(i, j)];
						result_same_imag += minus_pf_same_real * sigma_in[sigma_imag(i, j)]
						                  + minus_pf_same_imag * sigma_in[sigma_real(i, j)];
					}

					// bath b for site j
					const int b_j = BATHS_COUPLING[j * BATHS_MAX_PER_STATE + b];

					// check for coupling
					if (b_j >= 0) {
						// complex product: -pf_same[b_i] * sigma_in[id]
						const real_t minus_pf_same_real = -pf_same[2 * b_j];
						const real_t minus_pf_same_imag = -pf_same[2 * b_j + 1];
						result_same_real += minus_pf_same_real * sigma_in[sigma_real(i, j)]
						                  - minus_pf_same_imag * sigma_in[sigma_imag(i, j)];
						result_same_imag += minus_pf_same_real * sigma_in[sigma_imag(i, j)]
						                  + minus_pf_same_imag * sigma_in[sigma_real(i, j)];
					}
				} // for (b < BATHS_MAX_PER_STATE)
			} // if (i != j)
#endif // HEOM_ODE_DISABLE_SAME

			// -------------------- PLUS TERM --------------------
			real_t result_plus_real = 0.0;
			real_t result_plus_imag = 0.0;

#ifndef HEOM_ODE_DISABLE_PLUS
			// NOTE: code order of b_i and b_j plus terms changes the results
			//       so we use separate accumulation variables
			real_t result_plus_i_real = 0.0;
			real_t result_plus_i_imag = 0.0;
			real_t result_plus_j_real = 0.0;
			real_t result_plus_j_imag = 0.0;

			if (sigma_id < first_non_plus_id) {
				for (int b = 0; b < BATHS_MAX_PER_STATE; ++b) {
					// bath b for site i
					const int b_i = BATHS_COUPLING[i * BATHS_MAX_PER_STATE + b];

					// check for coupling
					if (b_i >= 0) {
						for (int m = 0; m < NUM_MATSUBARAS; ++m) {
							// Matsubara m for site i
							const int m_i = b_i * NUM_MATSUBARAS + m;

							// uid of ADO linked to site i
							// NOTE: can be accessed without checking ado_tuple_i below before due to
							//       first_non_plus_id check on top, this implies an order the host
							//       must guarantee
							const int plus_uid = ado_plus[sigma_id * ADO_WIDTH + m_i];
							if (mask[plus_uid]) {
								const real_t sigma_plus_site_i_real = sigma_in[sigma_uid_real(plus_uid, i, j)];
								const real_t sigma_plus_site_i_imag = sigma_in[sigma_uid_imag(plus_uid, i, j)];

								const int ado_tuple_i = ado_tuples[sigma_id * ADO_WIDTH + m_i];

								// NOTE: implicit complex i multiplication of scalar term
#ifdef HEOM_ODE_USE_SQRT_LUT
								const real_t aux_real = 0.0;
								const real_t aux_imag = SQRT_LUT[ado_tuple_i + 1] * cbk_sqrt[m_i] * CONST_SI_HBAR_INV_SQRT);
#else
								const real_t aux_real = 0.0;
								const real_t aux_imag = sqrt((ado_tuple_i + 1.0) * hypot(cbk[2 * m_i], cbk[2 * m_i + 1]) / CONST_SI_HBAR);
#endif
								// complex product of aux and sigma_plus_site_i_*
								// NOTE: full complex airthmetic including zero aux_real here for signed zero effects
								result_plus_i_real += aux_real * sigma_plus_site_i_real // NOTE: is zero
								                    - aux_imag * sigma_plus_site_i_imag;
								result_plus_i_imag += aux_real * sigma_plus_site_i_imag // NOTE: is zero
								                    + aux_imag * sigma_plus_site_i_real;
							} // if (mask[plus_uid])
						} // for (m < NUM_MATSUBARAS)
					} // if (b_i >= 0)

					// bath b for site j
					const int b_j = BATHS_COUPLING[j * BATHS_MAX_PER_STATE + b];

					// check for coupling
					if (b_j >= 0) {
						for (int m = 0; m < NUM_MATSUBARAS; ++m) {
							// Matsubara m for site j
							const int m_j = b_j * NUM_MATSUBARAS + m;

							// uid of ADO linked to site i
							// NOTE: can be accessed without checking ado_tuple_i below before due to
							//       first_non_plus_id check on top, this implies an order the host
							//       must guarantee
							const int plus_uid = ado_plus[sigma_id * ADO_WIDTH + m_j];
							if (mask[plus_uid]) {
								const real_t sigma_plus_site_j_real = sigma_in[sigma_uid_real(plus_uid, i, j)];
								const real_t sigma_plus_site_j_imag = sigma_in[sigma_uid_imag(plus_uid, i, j)];

								const int ado_tuple_j = ado_tuples[sigma_id * ADO_WIDTH + m_j];

								// (sqrt(..) / CONST_SI_HBAR) * COMPLEX_I
#ifdef HEOM_ODE_USE_SQRT_LUT
								const real_t aux_real = 0.0;
								const real_t aux_imag = SQRT_LUT[ado_tuple_j + 1] * cbk_sqrt[m_j] * CONST_SI_HBAR_INV_SQRT);
#else
								const real_t aux_real = 0.0;
								const real_t aux_imag = sqrt((ado_tuple_j + 1.0) * hypot(cbk[2 * m_j], cbk[2 * m_j + 1]) / CONST_SI_HBAR);
#endif
								// complex product of aux and sigma_plus_site_j_*
								// NOTE: full complex airthmetic including zero aux_real here for signed zero effects
								result_plus_j_real -= aux_real * sigma_plus_site_j_real // NOTE: is zero
								                    - aux_imag * sigma_plus_site_j_imag;
								result_plus_j_imag -= aux_real * sigma_plus_site_j_imag // NOTE: is zero
								                    + aux_imag * sigma_plus_site_j_real;
							} // if (mask[plus_uid])
						} // for (m < NUM_MATSUBARAS)
					} // if (b_j >= 0)

				} // for (b < BATHS_MAX_PER_STATE)
			} // if (sigma_id < first_non_plus_id)

			result_plus_real = result_plus_i_real + result_plus_j_real;
			result_plus_imag = result_plus_i_imag + result_plus_j_imag;
#endif // HEOM_ODE_DISABLE_PLUS

			// -------------------- MINUS TERM --------------------
			real_t result_minus_real = 0.0;
			real_t result_minus_imag = 0.0;

#ifndef HEOM_ODE_DISABLE_MINUS
			// NOTE: code order of b_i and b_j plus terms changes the results
			//       so we use separate accumulation variables
			real_t result_minus_i_real = 0.0;
			real_t result_minus_i_imag = 0.0;
			real_t result_minus_j_real = 0.0;
			real_t result_minus_j_imag = 0.0;

			for (int b = 0; b < BATHS_MAX_PER_STATE; ++b) {
				// bath b for site i
				const int b_i = BATHS_COUPLING[i * BATHS_MAX_PER_STATE + b];

				// check for coupling
				if (b_i >= 0) {
					for (int m = 0; m < NUM_MATSUBARAS; ++m) {
						// Matsubara m for site i
						const int m_i = b_i * NUM_MATSUBARAS + m;

						const int ado_tuple_i = ado_tuples[sigma_id * ADO_WIDTH + m_i];

						// NOTE: only access ado_minus with m_i if this condition holds
						if (ado_tuple_i > 0) {
							// uid of ADO linked to site i
							const int minus_uid = ado_minus[sigma_id * ADO_WIDTH + m_i];

							if (mask[minus_uid]) {
								const real_t sigma_minus_site_i_real = sigma_in[sigma_uid_real(minus_uid, i, j)];
								const real_t sigma_minus_site_i_imag = sigma_in[sigma_uid_imag(minus_uid, i, j)];

								const real_t cbk_m_i_real = cbk[2 * m_i];
								const real_t cbk_m_i_imag = cbk[2 * m_i + 1];

								// complex product: sigma_minus_site_i * cbk[m_i]
								const real_t prod_sigma_cbk_real = sigma_minus_site_i_real * cbk_m_i_real
								                                 - sigma_minus_site_i_imag * cbk_m_i_imag;
								const real_t prod_sigma_cbk_imag = sigma_minus_site_i_real * cbk_m_i_imag
								                                 + sigma_minus_site_i_imag * cbk_m_i_real;

								// complex product: complex i times product above, divided by CONST_SI_HBAR
								const real_t aux_real =  // 0.0 * prod_sigma_cbk_real
								                      - prod_sigma_cbk_imag  // -(1.0 * prod_sigma_cbk_imag)
								                      / CONST_SI_HBAR;
								const real_t aux_imag = // 0.0 * prod_sigma_cbk_imag
								                        prod_sigma_cbk_real  // (1.0 * prod_sigma_cbk_real)
								                      / CONST_SI_HBAR;

#ifdef HEOM_ODE_USE_SQRT_LUT
								const real_t sc = SQRT_LUT[ado_tuple_i] * cbk_sqrt[m_i] * CONST_SI_HBAR_SQRT;
#else
								const real_t sc = sqrt(ado_tuple_i / hypot(cbk_m_i_real, cbk_m_i_imag) * CONST_SI_HBAR);
#endif

								result_minus_i_real += aux_real * sc;
								result_minus_i_imag += aux_imag * sc;

								if (m == 0) {
									const real_t pf_circ_real = pf_circ[2 * b_i];
									const real_t pf_circ_imag = pf_circ[2 * b_i + 1];

									// complex product: (sigma_minus_site_i * pf_circ[b_i]) * sc
									result_minus_i_real += (  sigma_minus_site_i_real * pf_circ_real
									                        - sigma_minus_site_i_imag * pf_circ_imag
									                       ) * sc;
									result_minus_i_imag += (  sigma_minus_site_i_real * pf_circ_imag
									                        + sigma_minus_site_i_imag * pf_circ_real
									                       ) * sc;
								}
							} // if (mask[minus_uid])
						} // if (ado_tuple_i > 0)
					} // for (m < NUM_MATSUBARAS)
				} // if (b_i >= 0)

				// bath b for site j
				const int b_j = BATHS_COUPLING[j * BATHS_MAX_PER_STATE + b];

				// check for coupling
				if (b_j >= 0) {
					for (int m = 0; m < NUM_MATSUBARAS; ++m) {
						// Matsubara m for site i
						const int m_j = b_j * NUM_MATSUBARAS + m;

						const int ado_tuple_j = ado_tuples[sigma_id * ADO_WIDTH + m_j];
						
						// NOTE: only access ado_minus with m_j if this condition holds
						if (ado_tuple_j > 0) {
							// uid of ADO linked to site j
							const int minus_uid = ado_minus[sigma_id * ADO_WIDTH + m_j];

							if (mask[minus_uid]) {
								const real_t sigma_minus_site_j_real = sigma_in[sigma_uid_real(minus_uid, i, j)];
								const real_t sigma_minus_site_j_imag = sigma_in[sigma_uid_imag(minus_uid, i, j)];
								const real_t cbk_m_j_real = cbk[2 * m_j];
								const real_t cbk_m_j_imag = cbk[2 * m_j + 1];

								// complex product: sigma_minus_site_j * cbk[m_j]
								const real_t prod_sigma_cbk_real = sigma_minus_site_j_real * cbk_m_j_real
								                                 - sigma_minus_site_j_imag * cbk_m_j_imag;
								// NOTE: imaginary part seems a bit numerically unstable here (last 3-4 digits, maybe FMA-related)
								const real_t prod_sigma_cbk_imag = sigma_minus_site_j_real * cbk_m_j_imag
								                                 + sigma_minus_site_j_imag * cbk_m_j_real;

								const real_t aux_real =                      // 0.0 * prod_sigma_cbk_real
								                        prod_sigma_cbk_imag  // -(-1.0 * prod_sigma_cbk_imag)
								                      / CONST_SI_HBAR;
								const real_t aux_imag =                      // 0.0 * prod_sigma_cbk_imag
								                      - prod_sigma_cbk_real  // (-1.0 * prod_sigma_cbk_real)
								                      / CONST_SI_HBAR;

#ifdef HEOM_ODE_USE_SQRT_LUT
								const real_t sc = SQRT_LUT[ado_tuple_j] * cbk_sqrt[m_j] * CONST_SI_HBAR_SQRT;
#else
								const real_t sc = sqrt(ado_tuple_j / hypot(cbk_m_j_real, cbk_m_j_imag) * CONST_SI_HBAR);
#endif
								result_minus_j_real += aux_real * sc;
								result_minus_j_imag += aux_imag * sc;

								if (m == 0) {
									const real_t pf_circ_real = pf_circ[2 * b_j];
									const real_t pf_circ_imag = pf_circ[2 * b_j + 1];

									// complex product: (sigma_minus_site_j * pf_circ[b_j]) * sc
									result_minus_j_real += (  sigma_minus_site_j_real * pf_circ_real
									                        - sigma_minus_site_j_imag * pf_circ_imag
									                       ) * sc;
									result_minus_j_imag += (  sigma_minus_site_j_real * pf_circ_imag
									                        + sigma_minus_site_j_imag * pf_circ_real
									                       ) * sc;
								}
							} // if (mask[minus_uid])
						} // if (ado_tuple_j > 0)
					} // for (m < NUM_MATSUBARAS)
				} // if (b_j >= 0)

			} // for (b < BATHS_MAX_PER_STATE)

			result_minus_real = result_minus_i_real + result_minus_j_real;
			result_minus_imag = result_minus_i_imag + result_minus_j_imag;
#endif // HEOM_ODE_DISABLE_MINUS

			// --------------------OVERALL SUM OF TERMS --------------------

			real_t result_real = result_commutator_real + result_same_real + result_plus_real + result_minus_real;
			real_t result_imag = result_commutator_imag + result_same_imag + result_plus_imag + result_minus_imag;

#ifdef NOMA_NUM_ODE_ACCUMULATE
			if (acc_init) {
				// initialise
				sigma_acc[sigma_real(i, j)] = sigma_in[sigma_real(i, j)];
				sigma_acc[sigma_imag(i, j)] = sigma_in[sigma_imag(i, j)];
			}
			sigma_acc[sigma_real(i,j)] += acc_coeff * result_real; // do accumulation on the fly for y(n+1) = y(n) + h * b1 * k1 + ...
			sigma_acc[sigma_imag(i,j)] += acc_coeff * result_imag; // do accumulation on the fly for y(n+1) = y(n) + h * b1 * k1 + ...
#endif

#ifdef NOMA_NUM_ODE_SUBDIAGONAL
			// subdiagonal optimisation for writing the second argument for f' in: k2 = f'(..., y_n + h * a10 * k1), when computing k1
			sigma_acc[sigma_real(i,j)] = sigma_in[sigma_real(i,j)] + out_coeff * result_real;
			sigma_acc[sigma_imag(i,j)] = sigma_in[sigma_imag(i,j)] + out_coeff * result_imag;
#else
			// pure evaluation
			sigma_out[sigma_real(i,j)] = result_real;
			sigma_out[sigma_imag(i,j)] = result_imag;
#endif

#ifdef NOMA_NUM_ODE_TRACK_FLOWS
			flow_from_above[sigma_real(i,j)] += acc_coeff * result_plus_real;
			flow_from_above[sigma_imag(i,j)] += acc_coeff * result_plus_imag;
			flow_from_below[sigma_real(i,j)] += acc_coeff * result_minus_real;
			flow_from_below[sigma_imag(i,j)] += acc_coeff * result_minus_imag;
#endif

#if defined(HEOM_WORK_ITEM_GRANULARITY_MATRIX)
		} // for j
	} // for i
#elif defined(HEOM_WORK_ITEM_GRANULARITY_ELEMENT)
	// nothing to do
#else
// TODO: fail
#endif
}
