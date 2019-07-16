// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_population_dynamics_solver_hpp
#define heom_population_dynamics_solver_hpp

#include <iostream>
#include <functional>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <noma/memory/array_to_file.hpp>
#include <noma/num/make_stepper.hpp>
#include <noma/ocl/helper.hpp>

#include "heom/common.hpp"
#include "heom/config.hpp"
#include "heom/instance.hpp"
#include "heom/ocl_config.hpp"
#include "heom/rk_weighted_add.hpp"
#include "heom/solver_base.hpp"

namespace heom {

template<typename ODE_T, typename STEPPER_T, typename NORM_T>
class population_dynamics_solver : public solver_base<ODE_T, NORM_T>
{
public:
	population_dynamics_solver(const heom::ocl_config& ocl_config, const ocl::nd_range& range, const config& heom_config, instance& heom_instance);
	~population_dynamics_solver();

	// iterate
	void step_forward(size_t n);

	// statistics
	void output_debug_results(complex_t* buffer);
	statistics_vector stats();
	void write_ocl_runtime_config(std::ostream& os);
	void write_runtime_stats(std::ostream& os);
	void write_runtime_summary(std::ostream& os);

	// flow
	void update_flow_norm();
	void flow_norm(real_t* h_mem_flow_from_above_norm, real_t* h_mem_flow_from_below_norm);
	void reset_flows();

	// explicitly import super-class methods
	// avoids having to write 'this->' in front of them to enforce method lookup in dependant base class template
	using base_type = solver_base<ODE_T, NORM_T>;
	using base_type::heom_config;
	using base_type::heom_instance;
	using base_type::ocl_helper;
	using base_type::ode;
	using base_type::norm;
	using base_type::hierarchy_mask;
	using base_type::write_ocl_runtime_config;
	using base_type::runtime_stats_instance_values;
	using base_type::write_runtime_stats;
	using base_type::write_runtime_summary;
	using base_type::step_stats;
	using base_type::ocl_source_header;
	using base_type::current_hierarchy_buffer;
	using base_type::read_buffer;

	void reset();



	// NOTE: access needed for app_rate_kernel
	real_t time() const { return time_; }
	cl::Buffer& front_buffer() { return switch_buffers() ? d_mem_yn_b_ : d_mem_yn_a_; }
	cl::Buffer& back_buffer() { return switch_buffers() ? d_mem_yn_a_ : d_mem_yn_b_; }

private:
	void initialize();
	static const std::string ode_ocl_source_header(const config& heom_config); // write OpenCL ODE config to OpenCL source header

	bool switch_buffers() const	{ return step_count_ % 2 != 0; }

	// OpenCL kernels
	STEPPER_T stepper_;

	// OpenCL hierarchy buffers
	cl::Buffer d_mem_yn_a_;
	cl::Buffer d_mem_yn_b_;

	// OpenCL flow buffers
	complex_t* h_mem_flow_from_above_ = nullptr;
	complex_t* h_mem_flow_from_below_ = nullptr;
	cl::Buffer d_mem_flow_from_above_;
	cl::Buffer d_mem_flow_from_below_;
	cl::Buffer d_mem_flow_from_above_norm_;
	cl::Buffer d_mem_flow_from_below_norm_;

	// parameters (1-dim)
	real_t step_size_; // current step size
	std::vector<std::pair<real_t, size_t>> step_sizes_; // record of used step sizes as pair (step size, number of steps with that size), sum of first should equal time, sum of second should equal step_count, every change in step_size leads to a new record
	size_t step_count_; // overall step count
	real_t time_; // current time

};

template<typename ODE_T, typename STEPPER_T, typename NORM_T>
population_dynamics_solver<ODE_T, STEPPER_T, NORM_T>::population_dynamics_solver(const heom::ocl_config& ocl_config, const ocl::nd_range& range, const heom::config& heom_config, instance& heom_instance)
	: solver_base<ODE_T, NORM_T>(ocl_config, range, heom_config, heom_instance, ode_ocl_source_header(heom_config)),
	  stepper_(ocl_config.opencl_kernel_file_rk_weighted_add().empty() ?
	            num::make_stepper<ODE_T, STEPPER_T>(heom_config.solver_stepper_type(), ocl_helper(), heom::rk_weighted_add_ocl_source, heom::rk_weighted_add_ocl_kernel_name, ocl_source_header(), "", range, ode())
	          : num::make_stepper<ODE_T, STEPPER_T>(heom_config.solver_stepper_type(), ocl_helper(), boost::filesystem::path(ocl_config.opencl_kernel_file_rk_weighted_add()), ocl_config.opencl_kernel_name_rk_weighted_add(), ocl_source_header(), "", range, ode())),
	  step_size_(heom_config.solver_step_size()),
	  step_count_(0),
	  time_(0.0)
{
	initialize();
}

template<typename ODE_T, typename STEPPER_T, typename NORM_T>
population_dynamics_solver<ODE_T, STEPPER_T, NORM_T>::~population_dynamics_solver()
{
	if (heom_config().solver_track_flows()) {
		memory::aligned::instance().free(h_mem_flow_from_above_);
		memory::aligned::instance().free(h_mem_flow_from_below_);
	}
}

template<typename ODE_T, typename STEPPER_T, typename NORM_T>
void population_dynamics_solver<ODE_T, STEPPER_T, NORM_T>::initialize()
{
	step_sizes_.emplace_back(step_size_, 0);

	// allocate ode-state OpenCL memory
	d_mem_yn_a_ = ocl_helper().create_buffer(CL_MEM_READ_WRITE, heom_instance().size_hierarchy_byte(), nullptr);
	d_mem_yn_b_ = ocl_helper().create_buffer(CL_MEM_READ_WRITE, heom_instance().size_hierarchy_byte(), nullptr);

	// allocate flows memory, if configured
	if (heom_config().solver_track_flows()) {
		h_mem_flow_from_above_ = memory::aligned::instance().allocate<complex_t>(heom_instance().size_hierarchy_byte() / sizeof(complex_t));
		h_mem_flow_from_below_ = memory::aligned::instance().allocate<complex_t>(heom_instance().size_hierarchy_byte() / sizeof(complex_t));


		// allocate flow memory
		d_mem_flow_from_above_ = ocl_helper().create_buffer(CL_MEM_READ_WRITE, heom_instance().size_hierarchy_byte(), nullptr);
		d_mem_flow_from_below_ = ocl_helper().create_buffer(CL_MEM_READ_WRITE, heom_instance().size_hierarchy_byte(), nullptr);

		ode().set_flow_buffers(d_mem_flow_from_above_, d_mem_flow_from_below_);

		// allocate flow norm memory
		d_mem_flow_from_above_norm_ = ocl_helper().create_buffer(CL_MEM_READ_WRITE, heom_instance().matrices() * sizeof(real_t), nullptr);
		d_mem_flow_from_below_norm_ = ocl_helper().create_buffer(CL_MEM_READ_WRITE, heom_instance().matrices() * sizeof(real_t), nullptr);
	}

	reset();
}

template<typename ODE_T, typename STEPPER_T, typename NORM_T>
void population_dynamics_solver<ODE_T, STEPPER_T, NORM_T>::reset()
{
	cl_int err = 0;

	// write initial y_n
	err = ocl_helper().queue().enqueueWriteBuffer(d_mem_yn_a_, CL_TRUE, 0, heom_instance().size_hierarchy_byte(), heom_instance().hierarchy(), NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_yn_a_)");

	// re-init time-dependant members
	step_count_ = 0;
	time_ = 0.0;

	// NOTE: step_count is indirectly used here on the right side
	current_hierarchy_buffer() = front_buffer(); // set super class state
}

// result is passed to super class solver_base and then used for ODE_T construction
template<typename ODE_T, typename STEPPER_T, typename NORM_T>
const std::string population_dynamics_solver<ODE_T, STEPPER_T, NORM_T>::ode_ocl_source_header(const config& heom_config)
{
	std::stringstream options_stream;

	num::make_stepper_ode_compile_option<ODE_T, STEPPER_T>(heom_config.solver_stepper_type(), options_stream);

	if (heom_config.solver_track_flows())
	{
		options_stream << "#define HEOM_ODE_TRACK_FLOWS" << "\n";
	}

	return options_stream.str();
}

template<typename ODE_T, typename STEPPER_T, typename NORM_T>
void population_dynamics_solver<ODE_T, STEPPER_T, NORM_T>::update_flow_norm()
{
	// compute norm
	norm().update(d_mem_flow_from_above_norm_, d_mem_flow_from_above_);
	norm().update(d_mem_flow_from_below_norm_, d_mem_flow_from_below_);
}

template<typename ODE_T, typename STEPPER_T, typename NORM_T>
void population_dynamics_solver<ODE_T, STEPPER_T, NORM_T>::flow_norm(real_t* h_mem_flow_from_above_norm, real_t* h_mem_flow_from_below_norm)
{
	read_buffer(d_mem_flow_from_above_norm_, h_mem_flow_from_above_norm, heom_instance().matrices());
	read_buffer(d_mem_flow_from_below_norm_, h_mem_flow_from_below_norm, heom_instance().matrices());
}

template<typename ODE_T, typename STEPPER_T, typename NORM_T>
void population_dynamics_solver<ODE_T, STEPPER_T, NORM_T>::reset_flows()
{
	// set flow host memory to zero
	memset(h_mem_flow_from_above_, 0.0, sizeof(*h_mem_flow_from_above_));
	memset(h_mem_flow_from_below_, 0.0, sizeof(*h_mem_flow_from_below_));

	// write host to device memory
	cl_int err = 0;
	err = ocl_helper().queue().enqueueWriteBuffer(d_mem_flow_from_above_, CL_TRUE, 0, heom_instance().size_hierarchy_byte(), h_mem_flow_from_above_, NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_flow_from_above_)");
	err = ocl_helper().queue().enqueueWriteBuffer(d_mem_flow_from_below_, CL_TRUE, 0, heom_instance().size_hierarchy_byte(), h_mem_flow_from_below_, NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_flow_from_below_)");
}

template<typename ODE_T, typename STEPPER_T, typename NORM_T>
void population_dynamics_solver<ODE_T, STEPPER_T, NORM_T>::output_debug_results(complex_t* buffer)
{
	read_buffer(front_buffer(), buffer, heom_instance().size_hierarchy_top());
	write_complex_matrix(buffer, heom_instance().states(), std::cout);
	read_buffer(back_buffer(), buffer, heom_instance().size_hierarchy_top());
	write_complex_matrix(buffer, heom_instance().states(), std::cout);
}

template<typename ODE_T, typename STEPPER_T, typename NORM_T>
void population_dynamics_solver<ODE_T, STEPPER_T, NORM_T>::step_forward(size_t n)
{
	for (size_t i = 0; i < n; ++i) {
		bmt::timer t; // measure the whole step

		if (heom_config().solver_track_flows()) {
			reset_flows();
		}

		stepper_.step(time_, step_size_, front_buffer(), back_buffer()); // NOTE: front_buffer() is read, back_buffer() is written

		step_stats().add(t);

		// update time and step_size state and records
		++step_count_; // NOTE: front_buffer() and back_buffer() will switch
		current_hierarchy_buffer() = front_buffer(); // update super class state

		// TODO(adaptive time steps): get the actually used step size as a result from the stepper, and feed it back into it for the next iteration
		time_ += step_size_; // advance absolute time
		if( step_sizes_.back().first == step_size_ ) { // NOTE: intentional floating point equality here
			// step_size has not changed, increment counter for current record
			++(step_sizes_.back().second); // increment step_couter for that step_size
		} else {
			// new step_size, i.e. start a new record
			step_sizes_.emplace_back(step_size_, 1);
		}

		if (heom_config().solver_track_flows()) {
			// TODO: when writing the flow file, maybe sum up over multiple steps depending on observe_steps
			update_flow_norm();
		}
	}
}

template<typename ODE_T, typename STEPPER_T, typename NORM_T>
statistics_vector population_dynamics_solver<ODE_T, STEPPER_T, NORM_T>::stats()
{
	statistics_vector stats_vec;

	stats_vec.push_back(statistics_vector::value_type(step_stats(), "solver step"));
	stats_vec.push_back(statistics_vector::value_type(ode().kernel_stats(), "  heom_ode kernel"));
	stats_vec.push_back(statistics_vector::value_type(stepper_.kernel_stats(), "  rk_weighted_add kernel"));
	stats_vec.push_back(statistics_vector::value_type(norm().kernel_stats(), "  hierarchy_norm kernel"));

	return stats_vec;
}


template<typename ODE_T, typename STEPPER_T, typename NORM_T>
void population_dynamics_solver<ODE_T, STEPPER_T, NORM_T>::write_ocl_runtime_config(std::ostream& os)
{
	base_type::write_ocl_runtime_config(os);

	os << "OpenCL kernel for rk_weighted_add uses source from: " << (stepper_.uses_kernel_file() ? stepper_.kernel_file_name() : "[embedded]")
	   << std::endl;
}

template<typename ODE_T, typename STEPPER_T, typename NORM_T>
void population_dynamics_solver<ODE_T, STEPPER_T, NORM_T>::write_runtime_stats(std::ostream& os)
{
	base_type::write_runtime_stats(os);

	const std::string instance_values = runtime_stats_instance_values();

	// kernel_rk_weighted_add_stats
	os << "kernel_rk_weighted_add_stats" << default_delimiter
	   << stepper_.kernel_stats().string() << default_delimiter
	   << instance_values
	   << std::endl;
}

template<typename ODE_T, typename STEPPER_T, typename NORM_T>
void population_dynamics_solver<ODE_T, STEPPER_T, NORM_T>::write_runtime_summary(std::ostream& os)
{
	base_type::write_runtime_summary(os);

	write_time_statistics_summary(stepper_.kernel_stats(), "rk_weighted_add kernel", os);
}

} // namespace heom

#endif // heom_population_dynamics_solver_hpp
