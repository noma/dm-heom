// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_solver_base_hpp
#define heom_solver_base_hpp

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
#include "heom/hierarchy_mask.hpp"
#include "heom/instance.hpp"
#include "heom/ocl_config.hpp"
// TODO: remove, see below
#include "heom/hierarchy_arithmetic_mmult_left.hpp"
#include "heom/hierarchy_arithmetic_mmult_right.hpp"
#include "heom/hierarchy_arithmetic_scale.hpp"

namespace heom {

// TODO refactor inheriting classes to pass MMULT_* upwards, and remove includes and defeault arugments here
//template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
template<typename ODE_T, typename NORM_T, typename SCALE_T = heom::hierarchy_arithmetic_scale, typename MMULT_LEFT_T = heom::hierarchy_arithmetic_mmult_left, typename MMULT_RIGHT_T = heom::hierarchy_arithmetic_mmult_right>
class solver_base
{
public:
	solver_base(const heom::ocl_config& ocl_config, const ocl::nd_range& range, const config& heom_config, instance& heom_instance, const std::string& ode_ocl_source_header);
	~solver_base() = default;

	const heom::ocl_config& ocl_config() const { return ocl_config_; }
	const ocl::nd_range& range() const { return range_; }
	const config& heom_config() const { return heom_config_; }

	instance& heom_instance() { return heom_instance_; }
	ocl::helper& ocl_helper() { return ocl_helper_; }
	ODE_T& ode() { return ode_; }
	NORM_T& norm() { return norm_; }

	void update_instance();
	void update_hierarchy_norm();
	void update_hamiltonian(const real_t* new_hamiltonian);
	void update_hierarchy_top(const complex_t* buffer);

	// mask filtering
	void update_hierarchy_mask(const mask_t* new_mask_protection);
	void update_hierarchy_mask_protection(const mask_t* new_mask_protection);
	void update_hierarchy_mask_threshold(real_t threshold) { hierarchy_mask_threshold_ = threshold; }
	real_t hierarchy_mask_threshold() { return hierarchy_mask_threshold_; }
	cl::Buffer& hierarchy_mask() { return d_mem_hierarchy_mask_; }

	// data access
	void hierarchy(complex_t* buffer);
	void hierarchy_top(complex_t* buffer);
	void hierarchy_norm(real_t* buffer);
	void hierarchy_mask(mask_t* buffer);

	// general hierarchy arithmetic
	void hierarchy_mmult_left(const complex_t* matrix);
	void hierarchy_mmult_right(const complex_t* matrix);
	void hierarchy_scale(complex_t factor);
	void hierarchy_scale(cl::Buffer& d_mem_src, complex_t factor, cl::Buffer& d_mem_dst);

	// config and statistics
	void write_ocl_runtime_config(std::ostream& os);
	void write_runtime_stats(std::ostream& os);
	void write_runtime_summary(std::ostream& os);
	bmt::statistics& step_stats() { return step_stats_; }

	// buffer access
	cl::Buffer& current_hierarchy_buffer() { return d_mem_current_hierarchy_; }
	template<typename HOST_BUFFER_T>
	void read_buffer(cl::Buffer& d_mem_to_read, HOST_BUFFER_T* host_buffer, size_t size_byte, size_t offset = 0);

protected:
	const std::string ocl_source_header() const;  // generate OpenCL source header with defines needed by all kernels
	const std::string runtime_stats_instance_values();

private:
	void initialize();

	const heom::ocl_config& ocl_config_;
	const ocl::nd_range range_;
	const config& heom_config_;

	instance& heom_instance_;
	ocl::helper ocl_helper_;

	// OpenCL kernels
	ODE_T ode_;
	NORM_T norm_;
	SCALE_T scale_;
	MMULT_LEFT_T mmult_left_;
	MMULT_RIGHT_T mmult_right_;

	// OpenCL
	cl::Buffer d_mem_current_hierarchy_; // where the most recent hierarchy state can be found, to be updated by subclassses if needed

	// OpenCL norm and mask
	cl::Buffer d_mem_hierarchy_norm_;
	cl::Buffer d_mem_hierarchy_mask_;
	real_t hierarchy_mask_threshold_;

	// OpenCL temporary buffers for hierarchy arithmetic
	cl::Buffer d_mem_mmult_matrix_;

	// for benchmarking
	bmt::statistics step_stats_;
};

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::solver_base(const heom::ocl_config& ocl_config, const ocl::nd_range& range, const config& heom_config, instance& heom_instance, const std::string& ode_ocl_source_header)
	: ocl_config_(ocl_config), range_(range), heom_config_(heom_config), heom_instance_(heom_instance),
	  ocl_helper_(ocl_config),
	  ode_(ocl_config.opencl_kernel_file_heom_ode().empty() ? // was a kernel file name left unconfigured, i.e. empty
	        ODE_T(ocl_helper_, ocl_source_header() + ode_ocl_source_header, "", range, heom_instance_) // use ctor for embedded code
	      : ODE_T(ocl_helper_, ocl_config.opencl_kernel_file_heom_ode(), ocl_config.opencl_kernel_name_heom_ode(), ocl_source_header() + ode_ocl_source_header, "", range, heom_instance_) // provide configure path and kernel name
	     ),
	  norm_(ocl_config.opencl_kernel_file_hierarchy_norm().empty() ?
	         NORM_T(ocl_helper_, ocl_source_header(), "", range)
	       : NORM_T(ocl_helper_, ocl_config.opencl_kernel_file_hierarchy_norm(), ocl_config.opencl_kernel_name_hierarchy_norm(), ocl_source_header(), "", range)
	      ),
	  scale_(ocl_config.opencl_kernel_file_hierarchy_scale().empty() ?
	         SCALE_T(ocl_helper_, ocl_source_header(), "", range)
	       : SCALE_T(ocl_helper_, ocl_config.opencl_kernel_file_hierarchy_scale(), ocl_config.opencl_kernel_name_hierarchy_scale(), ocl_source_header(), "", range)
	  ),
	  mmult_left_(ocl_config.opencl_kernel_file_hierarchy_mmult_left().empty() ?
	              MMULT_LEFT_T(ocl_helper_, ocl_source_header(), "", range)
	            : MMULT_LEFT_T(ocl_helper_, ocl_config.opencl_kernel_file_hierarchy_mmult_left(), ocl_config.opencl_kernel_name_hierarchy_mmult_left(), ocl_source_header(), "", range)
	  ),
	  mmult_right_(ocl_config.opencl_kernel_file_hierarchy_mmult_right().empty() ?
	              MMULT_RIGHT_T(ocl_helper_, ocl_source_header(), "", range)
	            : MMULT_RIGHT_T(ocl_helper_, ocl_config.opencl_kernel_file_hierarchy_mmult_right(), ocl_config.opencl_kernel_name_hierarchy_mmult_right(), ocl_source_header(), "", range)
	  ),
	  step_stats_(heom_config.solver_steps())
{
	initialize();
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::initialize()
{
	cl_int err = 0;

	// allocate norm and mask OpenCL memory
	d_mem_hierarchy_norm_ = ocl_helper_.create_buffer(CL_MEM_READ_WRITE, heom_instance_.matrices() * sizeof(real_t), nullptr);
	d_mem_hierarchy_mask_ = ocl_helper_.create_buffer(CL_MEM_READ_WRITE, heom_instance_.matrices() * sizeof(mask_t), nullptr);

	d_mem_mmult_matrix_ = ocl_helper_.create_buffer(CL_MEM_READ_WRITE, heom_instance_.size_hierarchy_top_byte(), nullptr);

	// no threshold by setting negative value
	update_hierarchy_mask_threshold(-1.0);

	// write initial d_mem_mask
	mask_t* non_active_mask = memory::aligned::instance().allocate<mask_t>(heom_instance_.matrices());
	std::fill(non_active_mask, non_active_mask + heom_instance_.matrices(), mask_t::MASK_TRUE);

	err = ocl_helper_.queue().enqueueWriteBuffer(d_mem_hierarchy_mask_, CL_TRUE, 0, sizeof(mask_t) * heom_instance_.matrices(), non_active_mask, NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_hierarchy_mask_)");

	ode_.set_hierarchy_mask(d_mem_hierarchy_mask_); // NOTE: set parameter once, no need to update if memory is re-written with new mask

	memory::aligned::instance().free(non_active_mask);
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
const std::string solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::ocl_source_header() const
{
	std::stringstream baths_coupling_ss;
	// NOTE: sites to states conversion
	state_baths_coupling(heom_config_.baths_coupling(), heom_instance_.sites_to_states_mode()).print_flat(baths_coupling_ss);

	std::stringstream options_stream;
	options_stream
		<< "#define NUM_STATES " << heom_instance_.states() << "\n"
		<< "#define NUM_MATRICES " << heom_instance_.matrices() << "\n"
		<< "#define NUM_MATSUBARAS " << heom_instance_.matsubaras() << "\n"
		<< "#define BATHS_MAX_PER_STATE " << max_baths_per_state(heom_config_.baths_max_per_site(), heom_instance_.sites_to_states_mode()) << "\n"
		<< "#define BATHS_COUPLING_LIST " << baths_coupling_ss.str() << "\n"
		<< "#define ADO_WIDTH " << heom_instance_.ado_width() << "\n";

	return options_stream.str();
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
template<typename HOST_BUFFER_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::read_buffer(cl::Buffer& d_mem_to_read, HOST_BUFFER_T* host_buffer, size_t size, size_t offset)
{
	cl_int err = 0;
	err = ocl_helper_.queue().enqueueReadBuffer(d_mem_to_read, CL_TRUE, offset, size * sizeof(HOST_BUFFER_T), host_buffer, NULL, NULL);
	ocl::error_handler(err, "clEnqueueReadBuffer(d_mem_to_read)");
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::update_instance()
{
	read_buffer(d_mem_current_hierarchy_, heom_instance_.hierarchy(), heom_instance_.size_hierarchy_byte() / sizeof(real_t));
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::hierarchy(complex_t* buffer)
{
	read_buffer(d_mem_current_hierarchy_, buffer, heom_instance_.states() * heom_instance_.states() * heom_instance_.matrices());
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::hierarchy_top(complex_t* buffer)
{
	read_buffer(d_mem_current_hierarchy_, buffer, heom_instance_.states() * heom_instance_.states());
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::hierarchy_norm(real_t* buffer)
{
	read_buffer(d_mem_hierarchy_norm_, buffer, heom_instance_.matrices());
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::hierarchy_mask(mask_t* buffer)
{
	read_buffer(d_mem_hierarchy_mask_, buffer, heom_instance_.matrices());
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::hierarchy_scale(complex_t factor)
{
	scale_.run(d_mem_current_hierarchy_, factor, d_mem_current_hierarchy_);
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::hierarchy_scale(cl::Buffer& d_mem_src, complex_t factor, cl::Buffer& d_mem_dst)
{
	scale_.run(d_mem_src, factor, d_mem_dst);
}


template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::hierarchy_mmult_left(const complex_t* matrix)
{
	cl_int err = ocl_helper_.queue().enqueueWriteBuffer(d_mem_mmult_matrix_, CL_TRUE, 0,
	                                                    heom_instance_.size_hierarchy_top_byte(), matrix, NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_mmult_matrix_)");

	mmult_left_.run(d_mem_current_hierarchy_, d_mem_mmult_matrix_, d_mem_current_hierarchy_);
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::hierarchy_mmult_right(const complex_t* matrix)
{
	cl_int err = ocl_helper_.queue().enqueueWriteBuffer(d_mem_mmult_matrix_, CL_TRUE, 0,
	                                                    heom_instance_.size_hierarchy_top_byte(), matrix, NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_mmult_matrix_)");

	mmult_right_.run(d_mem_current_hierarchy_, d_mem_mmult_matrix_, d_mem_current_hierarchy_);
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::update_hierarchy_norm()
{
	norm_.update(hierarchy_mask_threshold_, d_mem_hierarchy_norm_, d_mem_hierarchy_mask_, d_mem_current_hierarchy_);
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::update_hamiltonian(const real_t* new_hamiltonian)
{
	heom_instance_.set_hamiltonian(new_hamiltonian);
	ode_.update_hamiltonian();
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::update_hierarchy_top(const complex_t* buffer)
{
	cl_int err = 0;
	err = ocl_helper_.queue().enqueueWriteBuffer(d_mem_current_hierarchy_, CL_TRUE, 0,
	                                             heom_instance_.size_hierarchy_top_byte(), buffer, NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_current_hierarchy_)");
}


template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::update_hierarchy_mask(const mask_t* new_mask)
{
	cl_int err = 0;
	err = ocl_helper_.queue().enqueueWriteBuffer(d_mem_hierarchy_mask_, CL_TRUE, 0,
	                                             sizeof(mask_t) * heom_instance_.matrices(), new_mask, NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_hierarchy_mask_)");
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::update_hierarchy_mask_protection(const mask_t* new_mask_protection)
{
	cl_int err = 0;

	mask_t* new_mask = memory::aligned::instance().allocate<mask_t>(heom_instance_.matrices());
	hierarchy_mask(new_mask); // write current mask into buffer

	for (int i = 0; i < heom_instance_.matrices(); i++) {
		if (new_mask_protection[i] == mask_t::MASK_PROTECTED)
			new_mask[i] = mask_t::MASK_PROTECTED;
	}

	err = ocl_helper_.queue().enqueueWriteBuffer(d_mem_hierarchy_mask_, CL_TRUE, 0,
	                                             sizeof(mask_t) * heom_instance_.matrices(), new_mask, NULL, NULL);
	ocl::error_handler(err, "clEnqueueWriteBuffer(d_mem_hierarchy_mask_)");
	memory::aligned::instance().free(new_mask);
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::write_ocl_runtime_config(std::ostream& os)
{
	os << "Using: " << std::endl;
	ocl_helper().write_device_info(os);
	os << std::endl;

	os << "Zero-copy OpenCL buffers configured:                " << (ocl_helper().zero_copy_available() ? "yes" : "no") << "\n"
	   << "OpenCL kernel for heom_ode uses source from:        " << (ode().uses_kernel_file() ? ode().kernel_file_name() : "[embedded]") << "\n"
	   << "OpenCL kernel for hierarchy_norm uses source from:  " << (norm().uses_kernel_file() ? norm().kernel_file_name() : "[embedded]")
	   << std::endl;
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
const std::string solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::runtime_stats_instance_values()
{
	std::stringstream instance_values;
	instance_values << heom_instance_.states() << default_delimiter
	                << heom_instance_.matrices() << default_delimiter
	                << heom_instance_.matsubaras() << default_delimiter
	                << heom_instance_.baths_per_state() << default_delimiter
	                << heom_instance_.ado_depth();

	return instance_values.str();
}


template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::write_runtime_stats(std::ostream& os)
{
	// header
	os << "name" << default_delimiter
	   << step_stats().header_string() << default_delimiter
	   << "states" << default_delimiter
	   << "matrices" << default_delimiter
	   << "matsubaras" << default_delimiter
	   << "baths_per_state" << default_delimiter
	   << "ado_depth" << '\n';

	const std::string instance_values = runtime_stats_instance_values();

	// step stats
	os << "solver step" << default_delimiter
	   << step_stats().string() << default_delimiter
	   << instance_values << '\n';

	// kernel_heom_ode_stats
	os << "kernel_heom_ode_stats" << default_delimiter
	   << ode().kernel_stats().string() << default_delimiter
	   << instance_values << '\n';

	// norm stats
	os << "kernel_hierarchy_norm_stats" << default_delimiter
	   << norm().kernel_stats().string() << default_delimiter
	   << instance_values << std::endl;
}

template<typename ODE_T, typename NORM_T, typename SCALE_T, typename MMULT_LEFT_T, typename MMULT_RIGHT_T>
void solver_base<ODE_T, NORM_T, SCALE_T, MMULT_LEFT_T, MMULT_RIGHT_T>::write_runtime_summary(std::ostream& os)
{
	write_time_statistics_summary(step_stats(), "solver step", os);
	write_time_statistics_summary(ode().kernel_stats(), "heom_ode kernel", os);
	write_time_statistics_summary(norm().kernel_stats(), "hierarchy_norm kernel", os);
}

} // namespace heom

#endif // heom_solver_base_hpp
