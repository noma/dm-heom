// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_thermal_state_search_solver_hpp
#define heom_thermal_state_search_solver_hpp

#include <iostream>
#include <functional>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <noma/memory/array_to_file.hpp>
#include <noma/ocl/helper.hpp>

#include "heom/common.hpp"
#include "heom/constants.hpp"
#include "heom/config.hpp"
#include "heom/eigen_decomposition.hpp"
#include "heom/instance.hpp"
#include "heom/ocl_config.hpp"
#include "heom/solver_base.hpp"

namespace heom {

template<typename ODE_T, typename TSS_T, typename NORM_T>
class thermal_state_search_solver : public solver_base<ODE_T, NORM_T>
{
public:
	thermal_state_search_solver(const heom::ocl_config& ocl_config, const ocl::nd_range& range, const config& heom_config, instance& heom_instance);
	~thermal_state_search_solver();

	// iterate
	void step_forward(size_t n);

	// check if we reached convergence
	bool test_convergence(real_t delta);

	// statistics
	void output_debug_results(complex_t* buffer);
	statistics_vector stats();
	void write_ocl_runtime_config(std::ostream& os);
	void write_runtime_stats(std::ostream& os);
	void write_runtime_summary(std::ostream& os);


	// explicitly import super-class methods
	// avoids having to write 'this->' in front of them to enforce method lookup in dependant base class template
	using base_type = solver_base<ODE_T, NORM_T>;
	using base_type::heom_config;
	using base_type::heom_instance;
	using base_type::ocl_helper;
	using base_type::ode;
	using base_type::norm;
	using base_type::update_hamiltonian;
	using base_type::hierarchy_mask;
	using base_type::write_ocl_runtime_config;
	using base_type::runtime_stats_instance_values;
	using base_type::write_runtime_stats;
	using base_type::write_runtime_summary;
	using base_type::step_stats;
	using base_type::ocl_source_header;
	using base_type::current_hierarchy_buffer;
	using base_type::read_buffer;
	using base_type::hierarchy_mmult_left;
	using base_type::hierarchy_mmult_right;

	void reset();

private:
	void initialize();
	static const std::string ode_ocl_source_header(const config& heom_config); // write OpenCL ODE config to OpenCL source header

	//bool switch_buffers() const	{ return step_count_ % 2 != 0; }
	bool switch_buffers() const	{ return false; } // NOTE: we never switch in this solver

	cl::Buffer& front_buffer() { return switch_buffers() ? d_mem_yn_b_ : d_mem_yn_a_; }
	cl::Buffer& back_buffer() { return switch_buffers() ? d_mem_yn_a_ : d_mem_yn_b_; }

	// OpenCL kernels
	TSS_T tss_;

	// OpenCL hierarchy buffers
	cl::Buffer d_mem_yn_a_;
	cl::Buffer d_mem_yn_b_;

	// parameters (1-dim)
	size_t step_count_; // overall step count

	// buffer for convergence check
	std::unique_ptr<complex_t[]> last_top_buffer_;
	std::unique_ptr<complex_t[]> current_top_buffer_;

	// diagonalising matrices of eigenvalue decomposition from formula (15)
	complex_matrix_t matrix_c_;
	complex_matrix_t matrix_c_transposed_;
};

template<typename ODE_T, typename TSS_T, typename NORM_T>
thermal_state_search_solver<ODE_T, TSS_T, NORM_T>::thermal_state_search_solver(const heom::ocl_config& ocl_config, const ocl::nd_range& range, const heom::config& heom_config, instance& heom_instance)
	: solver_base<ODE_T, NORM_T>(ocl_config, range, heom_config, heom_instance, ode_ocl_source_header(heom_config)),
	  tss_(ocl_config.opencl_kernel_file_thermal_state_search().empty() ?
	            TSS_T(ocl_helper(), ocl_source_header(), "", range)
	          : TSS_T(ocl_helper(), ocl_config.opencl_kernel_file_thermal_state_search(), ocl_config.opencl_kernel_name_thermal_state_search(), ocl_source_header(), "", range)
	  ),
	  step_count_(0),
	  last_top_buffer_(new complex_t[heom_instance.size_hierarchy_byte() / sizeof(complex_t)]),
	  current_top_buffer_(new complex_t[heom_instance.size_hierarchy_byte() / sizeof(complex_t)])
{
	initialize();
}

template<typename ODE_T, typename TSS_T, typename NORM_T>
thermal_state_search_solver<ODE_T, TSS_T, NORM_T>::~thermal_state_search_solver()
{
}

template<typename ODE_T, typename TSS_T, typename NORM_T>
void thermal_state_search_solver<ODE_T, TSS_T, NORM_T>::initialize()
{
	// setup thermal state search

	// expand hamiltonian if necessary
	auto hamiltonian = state_hamiltonian(heom_config().system_hamiltonian(), heom_instance().sites_to_states_mode());

	// prevent population of the ground state (if present) by assigning a high energy
	if (heom_instance().sites_to_states_mode() == sites_to_states_mode_t::with_ground_state) {
		hamiltonian.at(0, 0) = constants::tss_ground_state_energy;
	}
	DEBUG_ONLY( std::cout << "Hamiltonian: " << hamiltonian << std::endl; )

	// perform eigen decomposition of the Hamiltonian
	auto hamiltonian_eigen_decomp = heom::eigen_decomposition(hamiltonian);
	auto eigendiag = hamiltonian_eigen_decomp.eigendiag();
	auto eigenvectors = hamiltonian_eigen_decomp.eigenvectors();
	auto eigenvalues = hamiltonian_eigen_decomp.eigenvalues();

	DEBUG_ONLY( std::cout << "thermal_state_search_solver::initialize(..): eigenvalues: " << hamiltonian_eigen_decomp.eigenvalues() << std::endl; )
	DEBUG_ONLY( std::cout << "thermal_state_search_solver::initialize(..): eigenvectors: " << hamiltonian_eigen_decomp.eigenvectors() << std::endl; )
	DEBUG_ONLY( std::cout << "thermal_state_search_solver::initialize(..): diagonalised H: " << hamiltonian_eigen_decomp.eigendiag() << std::endl; )

	// Hamiltonian == eigenvectors * eigendiag * eigendiags^-1
	// eigendiag == hamiltonian_diag
	// hamiltonian_diag == eigenvectors^-1 * hamiltonian * eigenvectors

	// formula (15), (16): H_diag = C H C^T
	matrix_c_ = eigenvectors.transposed();
	matrix_c_transposed_ = eigenvectors;

	const real_t beta = 1.0 / (heom::constants::k_b * heom_config().baths_temperature());
	DEBUG_ONLY( std::cout << "thermal_state_search_solver::initialize(..): beta: " << beta << std::endl; )

	// diagnoalised hamiltonian H
	auto hamiltonian_diag = eigendiag;

	auto eigenvectors_map = eigen_decomposition::heom_matrix_to_eigen(eigenvectors);
	auto eigenvectors_transposed_map = eigen_decomposition::heom_matrix_to_eigen(eigenvectors.transposed());

	auto hamiltonian_map = eigen_decomposition::heom_matrix_to_eigen(hamiltonian);
	auto hamiltonian_diag_map = eigen_decomposition::heom_matrix_to_eigen(hamiltonian_diag);
	auto matrix_c_map = eigen_decomposition::heom_matrix_to_eigen(matrix_c_);
	auto matrix_c_transposed_map = eigen_decomposition::heom_matrix_to_eigen(matrix_c_transposed_);

	// Hamiltonian debug output
	DEBUG_ONLY( std::cout << "H hamiltonian:\n" << hamiltonian_map << '\n' << std::endl; )
	// Hamiltonian diagonalised
	DEBUG_ONLY( std::cout << "H^~ hamiltonian_diag:\n" << hamiltonian_diag_map << '\n' << std::endl; )
	// eigen test: use eigenvectors to diagonalise H
	DEBUG_ONLY( auto test_eigen = eigenvectors_transposed_map * hamiltonian_map * eigenvectors_map; )
	DEBUG_ONLY( std::cout << "H^~ = V^T * H * V:\n" << test_eigen << '\n' << std::endl; )
	// use formula (15) notation, should be the same
	DEBUG_ONLY( auto test_15 = matrix_c_map * hamiltonian_map * matrix_c_transposed_map; )
	DEBUG_ONLY( std::cout << "H^~ = C * H * C^T:\n" << test_15 << '\n' << std::endl; )

	// exponential form of H diag
	auto hamiltonian_diag_exp = hamiltonian_diag;
	complex_t diag_exp_sum = 0.0;
	for (size_t i = 0; i < hamiltonian_diag_exp.rows(); ++i) {
		hamiltonian_diag_exp.at(i, i) = std::exp(-hamiltonian_diag.at(i, i) * beta);
		diag_exp_sum += hamiltonian_diag_exp.at(i, i);
	}
	hamiltonian_diag_exp.scale(1.0/diag_exp_sum);

	auto hamiltonian_diag_exp_map = heom::eigen_decomposition::heom_matrix_to_eigen(hamiltonian_diag_exp);

	auto sigma_zero = matrix_c_transposed_map * hamiltonian_diag_exp_map * matrix_c_map;
	auto hierarchy_top = heom::eigen_decomposition::eigen_matrix_to_heom(sigma_zero);

	// hierarchy top
	DEBUG_ONLY( std::cout << "(16) sigma_0 = C^T * 'exp(H)' * C:\n"  << sigma_zero << '\n' << std::endl; )

	// initialise top matrix
	heom_instance().set_hierarchy_top(hierarchy_top.data());

	// setup pseudo hamiltonian with epsilon in every element
	// Version A: eps =  |E_N - E_1| / h_bar
	//const complex_t epsilon = std::abs(eigendiag.at(eigendiag.rows() - 1, eigendiag.cols() - 1) - eigendiag.at(0,0)) / heom::constants::h_bar;
	// Version B: eps = k_B * T
	//const complex_t epsilon = heom::constants::k_b * heom_config().baths_temperature();
	// Version C: magic number
	const complex_t epsilon = 1.0e14;
	DEBUG_ONLY( std::cout << "epsilon: " << epsilon << std::endl; )

	heom::complex_matrix_t pseudo_hamiltonian { hamiltonian_diag.rows(), hamiltonian_diag.cols(), epsilon };

	heom_instance().set_hamiltonian(reinterpret_cast<real_t*>(pseudo_hamiltonian.data()));
	update_hamiltonian(heom_instance().hamiltonian());

	// allocate ode-state OpenCL memory
	d_mem_yn_a_ = ocl_helper().create_buffer(CL_MEM_READ_WRITE, heom_instance().size_hierarchy_byte(), nullptr);
	d_mem_yn_b_ = ocl_helper().create_buffer(CL_MEM_READ_WRITE, heom_instance().size_hierarchy_byte(), nullptr);

	// setup tss kernel arguments
	tss_.set_args(heom_instance(), epsilon, matrix_c_, matrix_c_transposed_, eigenvalues);

	reset();
}

template<typename ODE_T, typename TSS_T, typename NORM_T>
void thermal_state_search_solver<ODE_T, TSS_T, NORM_T>::reset()
{
	cl_int err = 0;

	// write initial y_n
	err = ocl_helper().queue().enqueueWriteBuffer(d_mem_yn_a_, CL_TRUE, 0, heom_instance().size_hierarchy_byte(), heom_instance().hierarchy(), NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_yn_a_)");

	// re-init time-dependant members
	step_count_ = 0;

	// NOTE: step_count is indirectly used here on the right side
	current_hierarchy_buffer() = front_buffer(); // set super class state
}

// result is passed to super class solver_base and then used for ODE_T construction
template<typename ODE_T, typename TSS_T, typename NORM_T>
const std::string thermal_state_search_solver<ODE_T, TSS_T, NORM_T>::ode_ocl_source_header(const config& heom_config)
{
	std::stringstream options_stream;

	// configure thermal state search version of ODE kernel
	options_stream << "#define HEOM_TSS_EQUATION" << "\n";

	return options_stream.str();
}

template<typename ODE_T, typename TSS_T, typename NORM_T>
void thermal_state_search_solver<ODE_T, TSS_T, NORM_T>::output_debug_results(complex_t* buffer)
{
	read_buffer(front_buffer(), buffer, heom_instance().size_hierarchy_top_byte() / sizeof(complex_t));
	write_complex_matrix(buffer, heom_instance().states(), std::cout);
	read_buffer(back_buffer(), buffer, heom_instance().size_hierarchy_top_byte() / sizeof(complex_t));
	write_complex_matrix(buffer, heom_instance().states(), std::cout);
}

template<typename ODE_T, typename TSS_T, typename NORM_T>
void thermal_state_search_solver<ODE_T, TSS_T, NORM_T>::step_forward(size_t n)
{
	// copy last hierarchy top
	read_buffer(front_buffer(), last_top_buffer_.get(), heom_instance().size_hierarchy_top_byte() / sizeof(complex_t));

	for (size_t i = 0; i < n; ++i) {
		bmt::timer t; // measure the whole step

		// NOTE: front_buffer() and back_buffer() are manually alternated throughout one step here

		// call ODE formula (17) to (20)
		constexpr real_t step_size = 1.0;
		ode().solve(step_count_, step_size, front_buffer(), back_buffer());

		// every hierarchy node is now k_u, i.e. back_buffer()
		// formula (21) to (22)
		// (21) upper part
		current_hierarchy_buffer() = back_buffer(); // make sure we operate on the ODE result

		hierarchy_mmult_left(matrix_c_.data());
		hierarchy_mmult_right(matrix_c_transposed_.data());

		// (21) lower part
		tss_.run(back_buffer(), front_buffer()); // read from back_buffer, write into front_buffer

		// (22)
		current_hierarchy_buffer() = front_buffer(); // make sure we operate on the TSS result

		hierarchy_mmult_left(matrix_c_transposed_.data());
		hierarchy_mmult_right(matrix_c_.data());

		step_stats().add(t);

		// NOTE: result is in front buffer

		// update time and step_size state and records
		++step_count_; // NOTE: front_buffer() and back_buffer() do NOT switch implicitly here
	}
}

template<typename ODE_T, typename TSS_T, typename NORM_T>
bool thermal_state_search_solver<ODE_T, TSS_T, NORM_T>::test_convergence(real_t delta)
{
	// read current top buffer
	read_buffer(front_buffer(), current_top_buffer_.get(), heom_instance().size_hierarchy_top_byte() / sizeof(complex_t));

	// create norms as sum of absolute value of elements
	auto top_norm = [&](complex_t* top_buffer) {
		real_t norm = 0.0;
		const auto states = heom_instance().states();
		for (int_t i = 0; i < states; ++i) {
			for (int_t j = 0; j < states; ++j) {
				norm += std::abs(top_buffer[i * states + j]);
			}
		}
		return norm;
	};

	real_t last_norm = top_norm(last_top_buffer_.get());
	real_t current_norm = top_norm(current_top_buffer_.get());

	real_t diff = std::abs(last_norm - current_norm);

	DEBUG_ONLY( std::cout << "thermal_state_search_solver::test_convergence(delta = " << delta << "): diff = " << diff << std::endl; );
	return (diff < delta); // convergence criteria
}

template<typename ODE_T, typename TSS_T, typename NORM_T>
statistics_vector thermal_state_search_solver<ODE_T, TSS_T, NORM_T>::stats()
{
	statistics_vector stats_vec;

	stats_vec.push_back(statistics_vector::value_type(step_stats(), "solver step"));
	stats_vec.push_back(statistics_vector::value_type(ode().kernel_stats(), "  heom_ode kernel"));
	stats_vec.push_back(statistics_vector::value_type(tss_.kernel_stats(), "  thermal_state_search kernel"));
	stats_vec.push_back(statistics_vector::value_type(norm().kernel_stats(), "  hierarchy_norm kernel"));

	return stats_vec;
}


template<typename ODE_T, typename TSS_T, typename NORM_T>
void thermal_state_search_solver<ODE_T, TSS_T, NORM_T>::write_ocl_runtime_config(std::ostream& os)
{
	base_type::write_ocl_runtime_config(os);

	os << "OpenCL kernel for thermal_state_search uses source from: " << (tss_.uses_kernel_file() ? tss_.kernel_file_name() : "[embedded]")
	   << std::endl;
}

template<typename ODE_T, typename TSS_T, typename NORM_T>
void thermal_state_search_solver<ODE_T, TSS_T, NORM_T>::write_runtime_stats(std::ostream& os)
{
	base_type::write_runtime_stats(os);

	const std::string instance_values = runtime_stats_instance_values();

	// thermal_state_search stats
	os << "thermal_state_search_stats" << default_delimiter
	   << tss_.kernel_stats().string() << default_delimiter
	   << instance_values
	   << std::endl;
}

template<typename ODE_T, typename TSS_T, typename NORM_T>
void thermal_state_search_solver<ODE_T, TSS_T, NORM_T>::write_runtime_summary(std::ostream& os)
{
	base_type::write_runtime_summary(os);

	write_time_statistics_summary(tss_.kernel_stats(), "thermal_state_search kernel", os);
}

} // namespace heom

#endif // heom_thermal_state_search_solver_hpp
