// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_ode_hpp
#define heom_ode_hpp

#include <functional>

#include <noma/ocl/kernel_wrapper.hpp>

#include "heom/instance.hpp"

namespace heom {

class ode : public ocl::kernel_wrapper
{
public:
	// set of kernel arguments that usually do not change
	/*
	struct static_arg_set {
		std::pair<cl_uint, cl::Buffer> d_mem_hamiltonian;
		cl::Buffer d_mem_ado_tuples;
		cl::Buffer d_mem_ado_plus;
		cl::Buffer d_mem_ado_minus;
		cl::Buffer d_mem_pf_same_sum;
		cl::Buffer d_mem_pf_same;
		cl::Buffer d_mem_pf_circ;
		cl::Buffer d_mem_cbk;
		cl::Buffer d_med_mem_cbk_sqrt;
		int first_non_plus_id;
	};*/

	ode(ocl::helper& ocl, const std::string& source_header, const std::string& ocl_compile_options, const ocl::nd_range& range, instance& heom_instance);
	ode(ocl::helper& ocl, const boost::filesystem::path& heom_ode_file_name, const std::string& heom_ode_kernel_name, const std::string& source_header, const std::string& ocl_compile_options, const ocl::nd_range& range, instance& heom_instance);

	void solve(real_t time, real_t step_size, cl::Buffer& in, cl::Buffer& out);
	void solve(real_t time, real_t step_size, cl::Buffer& in, cl::Buffer& out, real_t acc_coeff);
	void solve(real_t time, real_t step_size, cl::Buffer& in, cl::Buffer& out, cl::Buffer& acc, real_t acc_coff, bool acc_init);
	void solve(real_t time, real_t step_size, cl::Buffer& in, cl::Buffer& out, real_t out_coeff, cl::Buffer& acc, real_t acc_coff, bool acc_init);

	size_t buffer_size_byte() { return heom_instance.size_hierarchy_byte(); }

	// get assigned heom instance from the ODE
	instance& get_heom_instance() { return heom_instance; }

	// re-transfer hamiltonian from instance class (assuming it was updated)
	void update_hamiltonian();

	void set_hierarchy_mask(cl::Buffer& d_mem_hierarchy_mask);

	void set_flow_buffers(cl::Buffer& flow_from_above, cl::Buffer& flow_from_below);

	using pre_evaluate_action_t = std::function<void(ode&, real_t, real_t, cl::Buffer&, cl::Buffer&)>;
	void add_pre_evaluate_action(const pre_evaluate_action_t& action) { pre_evaluate_actions.emplace_back(action); }
	void clear_pre_evaluate_actions() { pre_evaluate_actions.clear(); };

private:
	void initialise();

	void perform_pre_evaluate_actions(real_t time, real_t step_size, cl::Buffer& in, cl::Buffer& out);

	instance& heom_instance;

	// list of function to execute before evaluating the ODE, e.g. compute some time-dependant input
	std::vector<pre_evaluate_action_t> pre_evaluate_actions;

	// OpenCL buffers
	cl::Buffer d_mem_hamiltonian;
	cl::Buffer d_mem_ado_tuples;
	cl::Buffer d_mem_ado_plus;
	cl::Buffer d_mem_ado_minus;
	cl::Buffer d_mem_pf_same_sum;
	cl::Buffer d_mem_pf_same;
	cl::Buffer d_mem_pf_circ;
	cl::Buffer d_mem_cbk;
	cl::Buffer d_mem_cbk_sqrt;
	static const std::string embedded_ocl_source_;
	static const std::string embedded_ocl_kernel_name_;
};

} // namespace heom

#endif // heom_ode_hpp
