// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_ocl_config_hpp
#define heom_ocl_config_hpp

#include <string>

#include <boost/program_options.hpp>
#include <CL/cl.hpp>
#include <noma/ocl/config.hpp>

#include "heom/common.hpp"

namespace heom {

class ocl_config : public ocl::config
{
public:
	ocl_config(const std::string& config_file_name, bool parse_config = true);

	const std::string& opencl_kernel_file_heom_ode() const
	{ return opencl_kernel_file_heom_ode_; }

	const std::string& opencl_kernel_name_heom_ode() const
	{ return opencl_kernel_name_heom_ode_; }

	const std::string& opencl_kernel_file_rk_weighted_add() const
	{ return opencl_kernel_file_rk_weighted_add_; }

	const std::string& opencl_kernel_name_rk_weighted_add() const
	{ return opencl_kernel_name_rk_weighted_add_; }

	const std::string& opencl_kernel_file_thermal_state_search() const
	{ return opencl_kernel_file_thermal_state_search_; }

	const std::string& opencl_kernel_name_thermal_state_search() const
	{ return opencl_kernel_name_thermal_state_search_; }

	const std::string& opencl_kernel_file_hierarchy_norm() const
	{ return opencl_kernel_file_hierarchy_norm_; }

	const std::string& opencl_kernel_name_hierarchy_norm() const
	{ return opencl_kernel_name_hierarchy_norm_; }

	const std::string& opencl_kernel_file_error_norm() const
	{ return opencl_kernel_file_error_norm_; }

	const std::string& opencl_kernel_name_error_norm() const
	{ return opencl_kernel_name_error_norm_; }

	const std::string& opencl_kernel_file_hierarchy_add() const
	{ return opencl_kernel_file_hierarchy_add_; }

	const std::string& opencl_kernel_name_hierarchy_add() const
	{ return opencl_kernel_name_hierarchy_add_; }

	const std::string& opencl_kernel_file_hierarchy_subtract() const
	{ return opencl_kernel_file_hierarchy_subtract_; }

	const std::string& opencl_kernel_name_hierarchy_subtract() const
	{ return opencl_kernel_name_hierarchy_subtract_; }

	const std::string& opencl_kernel_file_hierarchy_scale() const
	{ return opencl_kernel_file_hierarchy_scale_; }

	const std::string& opencl_kernel_name_hierarchy_scale() const
	{ return opencl_kernel_name_hierarchy_scale_; }

	const std::string& opencl_kernel_file_hierarchy_dot_product() const
	{ return opencl_kernel_file_hierarchy_dot_product_; }

	const std::string& opencl_kernel_name_hierarchy_dot_product() const
	{ return opencl_kernel_name_hierarchy_dot_product_; }

	const std::string& opencl_kernel_file_hierarchy_mmult_left() const
	{ return opencl_kernel_file_hierarchy_mmult_left_; }

	const std::string& opencl_kernel_name_hierarchy_mmult_left() const
	{ return opencl_kernel_name_hierarchy_mmult_left_; }

	const std::string& opencl_kernel_file_hierarchy_mmult_right() const
	{ return opencl_kernel_file_hierarchy_mmult_right_; }

	const std::string& opencl_kernel_name_hierarchy_mmult_right() const
	{ return opencl_kernel_name_hierarchy_mmult_right_; }

private:
	std::string opencl_kernel_file_heom_ode_;
	std::string opencl_kernel_name_heom_ode_;
	std::string opencl_kernel_file_rk_weighted_add_;
	std::string opencl_kernel_name_rk_weighted_add_;
	std::string opencl_kernel_file_thermal_state_search_;
	std::string opencl_kernel_name_thermal_state_search_;
	std::string opencl_kernel_file_hierarchy_norm_;
	std::string opencl_kernel_name_hierarchy_norm_;
	std::string opencl_kernel_file_error_norm_;
	std::string opencl_kernel_name_error_norm_;
	std::string opencl_kernel_file_hierarchy_add_;
	std::string opencl_kernel_name_hierarchy_add_;
	std::string opencl_kernel_file_hierarchy_subtract_;
	std::string opencl_kernel_name_hierarchy_subtract_;
	std::string opencl_kernel_file_hierarchy_scale_;
	std::string opencl_kernel_name_hierarchy_scale_;
	std::string opencl_kernel_file_hierarchy_dot_product_;
	std::string opencl_kernel_name_hierarchy_dot_product_;
	std::string opencl_kernel_file_hierarchy_mmult_left_;
	std::string opencl_kernel_name_hierarchy_mmult_left_;
	std::string opencl_kernel_file_hierarchy_mmult_right_;
	std::string opencl_kernel_name_hierarchy_mmult_right_;
};

} // namespace heom

#endif // heom_ocl_config_hpp
