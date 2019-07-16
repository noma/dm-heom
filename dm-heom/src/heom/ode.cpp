// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/ode.hpp"

namespace heom {

const std::string ode::embedded_ocl_source_ {
#include "heom_ode.cl.hpp"  // NOTE: generated by CMake
};
const std::string ode::embedded_ocl_kernel_name_ { "heom_ode" };

void ode::perform_pre_evaluate_actions(real_t time, real_t step_size, cl::Buffer& in, cl::Buffer& out)
{
	for (auto& a : pre_evaluate_actions)
		a(*this, time, step_size, in, out); // call action
}

ode::ode(ocl::helper& ocl, const std::string& source_header, const std::string& ocl_compile_options, const ocl::nd_range& range, instance& heom_instance)
	: ocl::kernel_wrapper(ocl, embedded_ocl_source_, embedded_ocl_kernel_name_, source_header, ocl_compile_options, range), heom_instance(heom_instance)
{
	initialise();
}

ode::ode(ocl::helper& ocl, const boost::filesystem::path& heom_ode_file_name, const std::string& heom_ode_kernel_name, const std::string& source_header, const std::string& ocl_compile_options, const ocl::nd_range& range, instance& heom_instance)
	: ocl::kernel_wrapper(ocl, heom_ode_file_name, heom_ode_kernel_name, source_header, ocl_compile_options, range), heom_instance(heom_instance)
{
	initialise();
}

void ode::initialise()
{
	cl_int err = 0;

	// allocate OpenCL buffers
	d_mem_hamiltonian = ocl_.create_buffer(CL_MEM_READ_ONLY, heom_instance.size_hamiltonian_byte(), nullptr);
	d_mem_ado_tuples = ocl_.create_buffer(CL_MEM_READ_ONLY, heom_instance.size_ado_tuples_byte(), nullptr);
	d_mem_ado_plus = ocl_.create_buffer(CL_MEM_READ_ONLY, heom_instance.size_ado_plus_byte(), nullptr);
	d_mem_ado_minus = ocl_.create_buffer(CL_MEM_READ_ONLY, heom_instance.size_ado_minus_byte(), nullptr);
	d_mem_pf_same_sum  = ocl_.create_buffer(CL_MEM_READ_ONLY, heom_instance.size_pf_same_sum_byte(), nullptr);
	d_mem_pf_same = ocl_.create_buffer(CL_MEM_READ_ONLY, heom_instance.size_pf_same_byte(), nullptr);
	d_mem_pf_circ = ocl_.create_buffer(CL_MEM_READ_ONLY, heom_instance.size_pf_circ_byte(), nullptr);
	d_mem_cbk = ocl_.create_buffer(CL_MEM_READ_ONLY, heom_instance.size_cbk_byte(), nullptr);
	d_mem_cbk_sqrt = ocl_.create_buffer(CL_MEM_READ_ONLY, heom_instance.size_cbk_sqrt_byte(), nullptr);

	// inititialise OpenCL buffers
	err = ocl_.queue().enqueueWriteBuffer(d_mem_hamiltonian, CL_TRUE, 0, heom_instance.size_hamiltonian_byte(), heom_instance.hamiltonian(), NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_hamiltonian)");
	err = ocl_.queue().enqueueWriteBuffer(d_mem_ado_tuples, CL_TRUE, 0, heom_instance.size_ado_tuples_byte(), heom_instance.ado_tuples(), NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_ado_tuples)");
	err = ocl_.queue().enqueueWriteBuffer(d_mem_ado_plus, CL_TRUE, 0, heom_instance.size_ado_plus_byte(), heom_instance.ado_plus(), NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_ado_plus)");
	err = ocl_.queue().enqueueWriteBuffer(d_mem_ado_minus, CL_TRUE, 0, heom_instance.size_ado_minus_byte(), heom_instance.ado_minus(), NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_ado_minus)");
	err = ocl_.queue().enqueueWriteBuffer(d_mem_pf_same_sum, CL_TRUE, 0, heom_instance.size_pf_same_sum_byte(), heom_instance.pf_same_sum(), NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_pf_same_sum)");
	err = ocl_.queue().enqueueWriteBuffer(d_mem_pf_same, CL_TRUE, 0, heom_instance.size_pf_same_byte(), heom_instance.pf_same(), NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_pf_same)");
	err = ocl_.queue().enqueueWriteBuffer(d_mem_pf_circ, CL_TRUE, 0, heom_instance.size_pf_circ_byte(), heom_instance.pf_circ(), NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_pf_circ)");
	err = ocl_.queue().enqueueWriteBuffer(d_mem_cbk, CL_TRUE, 0, heom_instance.size_cbk_byte(), heom_instance.cbk(), NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_cbk)");
	err = ocl_.queue().enqueueWriteBuffer(d_mem_cbk_sqrt, CL_TRUE, 0, heom_instance.size_cbk_sqrt_byte(), heom_instance.cbk_sqrt(), NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_cbk_sqrt)");

	// set static kernel arguments
	err = kernel_.setArg(7, d_mem_hamiltonian);
	ocl::error_handler(err, "clSetKernelArg(7)");
	err = kernel_.setArg(8, d_mem_ado_tuples);
	ocl::error_handler(err, "clSetKernelArg(8)");
	err = kernel_.setArg(9, d_mem_ado_plus);
	ocl::error_handler(err, "clSetKernelArg(9)");
	err = kernel_.setArg(10, d_mem_ado_minus);
	ocl::error_handler(err, "clSetKernelArg(10)");
	err = kernel_.setArg(11, d_mem_pf_same_sum);
	ocl::error_handler(err, "clSetKernelArg(11)");
	err = kernel_.setArg(12, d_mem_pf_same);
	ocl::error_handler(err, "clSetKernelArg(12)");
	err = kernel_.setArg(13, d_mem_pf_circ);
	ocl::error_handler(err, "clSetKernelArg(13)");
	err = kernel_.setArg(14, d_mem_cbk);
	ocl::error_handler(err, "clSetKernelArg(14)");
	err = kernel_.setArg(15, d_mem_cbk_sqrt);
	ocl::error_handler(err, "clSetKernelArg(15)");
	err = kernel_.setArg(16, heom_instance.first_non_plus_id());
	ocl::error_handler(err, "clSetKernelArg(16)");
	err = kernel_.setArg(18, nullptr);
	ocl::error_handler(err, "clSetKernelArg(18)");
	err = kernel_.setArg(19, nullptr);
	ocl::error_handler(err, "clSetKernelArg(19)");
}

void ode::solve(real_t time, real_t step_size, cl::Buffer& in, cl::Buffer& out)
{
	perform_pre_evaluate_actions(time, step_size, in, out);
	cl_int err = 0;
	//std::cout << "ODE: time = " << time << " | step_size = " << step_size << std::endl;
	// set dynamic kernel arguments
	err = kernel_.setArg(0, step_size);
	ocl::error_handler(err, "clSetKernelArg(0)");
	err = kernel_.setArg(1, out);
	ocl::error_handler(err, "clSetKernelArg(1)");
	err = kernel_.setArg(2, 1.0);
	ocl::error_handler(err, "clSetKernelArg(2)");
	err = kernel_.setArg(3, NULL);
	ocl::error_handler(err, "clSetKernelArg(3)");
	err = kernel_.setArg(4, 1.0);
	ocl::error_handler(err, "clSetKernelArg(4)");
	err = kernel_.setArg(5, false);
	ocl::error_handler(err, "clSetKernelArg(5)");
	err = kernel_.setArg(6, in);
	ocl::error_handler(err, "clSetKernelArg(6)");

	// run and benchmark kernel
	run_kernel();
}

void ode::solve(real_t time, real_t step_size, cl::Buffer& in, cl::Buffer& out, real_t acc_coeff)
{
	perform_pre_evaluate_actions(time, step_size, in, out);
	cl_int err = 0;
	//std::cout << "ODE: time = " << time << " | step_size = " << step_size << std::endl;
	// set dynamic kernel arguments
	err = kernel_.setArg(0, step_size);
	ocl::error_handler(err, "clSetKernelArg(0)");
	err = kernel_.setArg(1, out);
	ocl::error_handler(err, "clSetKernelArg(1)");
	err = kernel_.setArg(2, 1.0);
	ocl::error_handler(err, "clSetKernelArg(2)");
	err = kernel_.setArg(3, NULL);
	ocl::error_handler(err, "clSetKernelArg(3)");
	err = kernel_.setArg(4, acc_coeff);
	ocl::error_handler(err, "clSetKernelArg(4)");
	err = kernel_.setArg(5, false);
	ocl::error_handler(err, "clSetKernelArg(5)");
	err = kernel_.setArg(6, in);
	ocl::error_handler(err, "clSetKernelArg(6)");

	// run and benchmark kernel
	run_kernel();
}

void ode::solve(real_t time, real_t step_size, cl::Buffer& in, cl::Buffer& out, cl::Buffer& acc, real_t acc_coff, bool acc_init)
{
	perform_pre_evaluate_actions(time, step_size, in, out);
	cl_int err = 0;
	//std::cout << "ODE: time = " << time << " | step_size = " << step_size << std::endl;
	// set dynamic kernel arguments
	err = kernel_.setArg(0, step_size);
	ocl::error_handler(err, "clSetKernelArg(0)");
	err = kernel_.setArg(1, out);
	ocl::error_handler(err, "clSetKernelArg(1)");
	err = kernel_.setArg(2, 1.0);
	ocl::error_handler(err, "clSetKernelArg(2)");
	err = kernel_.setArg(3, acc);
	ocl::error_handler(err, "clSetKernelArg(3)");
	err = kernel_.setArg(4, acc_coff);
	ocl::error_handler(err, "clSetKernelArg(4)");
	err = kernel_.setArg(5, acc_init);
	ocl::error_handler(err, "clSetKernelArg(5)");
	err = kernel_.setArg(6, in);
	ocl::error_handler(err, "clSetKernelArg(6)");

	// run and benchmark kernel
	run_kernel();
}

void ode::solve(real_t time, real_t step_size, cl::Buffer& in, cl::Buffer& out, real_t out_coeff, cl::Buffer& acc, real_t acc_coff, bool acc_init)
{
	perform_pre_evaluate_actions(time, step_size, in, out);
	cl_int err = 0;
	//std::cout << "ODE: time = " << time << " | step_size = " << step_size << std::endl;
	// set dynamic kernel arguments
	err = kernel_.setArg(0, step_size);
	ocl::error_handler(err, "clSetKernelArg(0)");
	err = kernel_.setArg(1, out);
	ocl::error_handler(err, "clSetKernelArg(1)");
	err = kernel_.setArg(2, out_coeff);
	ocl::error_handler(err, "clSetKernelArg(2)");
	err = kernel_.setArg(3, acc);
	ocl::error_handler(err, "clSetKernelArg(3)");
	err = kernel_.setArg(4, acc_coff);
	ocl::error_handler(err, "clSetKernelArg(4)");
	err = kernel_.setArg(5, acc_init);
	ocl::error_handler(err, "clSetKernelArg(5)");
	err = kernel_.setArg(6, in);
	ocl::error_handler(err, "clSetKernelArg(6)");

	// run and benchmark kernel
	run_kernel();
}


void ode::update_hamiltonian()
{
	// write new buffer
	cl_int err = 0;
	err = ocl_.queue().enqueueWriteBuffer(d_mem_hamiltonian, CL_TRUE, 0, heom_instance.size_hamiltonian_byte(), heom_instance.hamiltonian(), NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_hamiltonian)");
}

void ode::set_hierarchy_mask(cl::Buffer& d_mem_hierarchy_mask)
{
	cl_int err = 0;
	err = kernel_.setArg(17, d_mem_hierarchy_mask);
	ocl::error_handler(err, "clSetKernelArg(17)");
}

void ode::set_flow_buffers(cl::Buffer& d_mem_flow_from_above, cl::Buffer& d_mem_flow_from_below)
{
	cl_int err = 0;
	err = kernel_.setArg(18, d_mem_flow_from_above);
	ocl::error_handler(err, "clSetKernelArg(18)");
	err = kernel_.setArg(19, d_mem_flow_from_below);
	ocl::error_handler(err, "clSetKernelArg(19)");
}

} // namespace heom
