// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/ocl_config.hpp"

#include <fstream>
#include <iostream>

namespace heom {

namespace bpo = ::boost::program_options;

ocl_config::ocl_config(const std::string& config_file_name, bool parse_config)
	: ocl::config(config_file_name, false)
{
	// initialise program options description
	desc_.add_options()
		("opencl.kernel_file_heom_ode", bpo::value(&opencl_kernel_file_heom_ode_)->default_value(opencl_kernel_file_heom_ode_), "Optional path to the heom_ode OpenCL kernel source code file.")
		("opencl.kernel_name_heom_ode", bpo::value(&opencl_kernel_name_heom_ode_)->default_value(opencl_kernel_name_heom_ode_), "Function name of the heom_ode OpenCL kernel.")
		("opencl.kernel_file_rk_weighted_add", bpo::value(&opencl_kernel_file_rk_weighted_add_)->default_value(opencl_kernel_file_rk_weighted_add_), "Optional path to the rk_weighted_add OpenCL kernel source code file.")
		("opencl.kernel_name_rk_weighted_add", bpo::value(&opencl_kernel_name_rk_weighted_add_)->default_value(opencl_kernel_name_rk_weighted_add_), "Function name of the rk_weighted_add OpenCL kernel.")
		("opencl.kernel_file_hierarchy_norm", bpo::value(&opencl_kernel_file_hierarchy_norm_)->default_value(opencl_kernel_file_hierarchy_norm_), "Optional path to the hierarchy_norm OpenCL kernel source code file.")
		("opencl.kernel_name_hierarchy_norm", bpo::value(&opencl_kernel_name_hierarchy_norm_)->default_value(opencl_kernel_name_hierarchy_norm_), "Function name of the hierarchy_norm OpenCL kernel.")
		("opencl.kernel_file_error_norm", bpo::value(&opencl_kernel_file_error_norm_)->default_value(opencl_kernel_file_error_norm_), "Optional path to the error_norm OpenCL kernel source code file.")
		("opencl.kernel_name_error_norm", bpo::value(&opencl_kernel_name_error_norm_)->default_value(opencl_kernel_name_error_norm_), "Function name of the error_norm OpenCL kernel.")
		("opencl.kernel_file_hierarchy_add", bpo::value(&opencl_kernel_file_hierarchy_add_)->default_value(opencl_kernel_file_hierarchy_add_), "Optional path to the hierarchy_add OpenCL kernel source code file.")
		("opencl.kernel_name_hierarchy_add", bpo::value(&opencl_kernel_name_hierarchy_add_)->default_value(opencl_kernel_name_hierarchy_add_), "Function name of the hierarchy_add OpenCL kernel.")
		("opencl.kernel_file_hierarchy_subtract", bpo::value(&opencl_kernel_file_hierarchy_subtract_)->default_value(opencl_kernel_file_hierarchy_subtract_), "Optional path to the hierarchy_subtract OpenCL kernel source code file.")
		("opencl.kernel_name_hierarchy_subtract", bpo::value(&opencl_kernel_name_hierarchy_subtract_)->default_value(opencl_kernel_name_hierarchy_subtract_), "Function name of the hierarchy_subtract OpenCL kernel.")
		("opencl.kernel_file_hierarchy_scale", bpo::value(&opencl_kernel_file_hierarchy_scale_)->default_value(opencl_kernel_file_hierarchy_scale_), "Optional path to the hierarchy_scale OpenCL kernel source code file.")
		("opencl.kernel_name_hierarchy_scale", bpo::value(&opencl_kernel_name_hierarchy_scale_)->default_value(opencl_kernel_name_hierarchy_scale_), "Function name of the hierarchy_scale OpenCL kernel.")
		("opencl.kernel_file_hierarchy_dot_product", bpo::value(&opencl_kernel_file_hierarchy_dot_product_)->default_value(opencl_kernel_file_hierarchy_dot_product_), "Optional path to the hierarchy_dot_product OpenCL kernel source code file.")
		("opencl.kernel_name_hierarchy_dot_product", bpo::value(&opencl_kernel_name_hierarchy_dot_product_)->default_value(opencl_kernel_name_hierarchy_dot_product_), "Function name of the hierarchy_dot_product OpenCL kernel.")
		("opencl.kernel_file_hierarchy_mmult_left", bpo::value(&opencl_kernel_file_hierarchy_mmult_left_)->default_value(opencl_kernel_file_hierarchy_mmult_left_), "Optional path to the hierarchy_mmult_left OpenCL kernel source code file.")
		("opencl.kernel_name_hierarchy_mmult_left", bpo::value(&opencl_kernel_name_hierarchy_mmult_left_)->default_value(opencl_kernel_name_hierarchy_mmult_left_), "Function name of the hierarchy_mmult_left OpenCL kernel.")
		("opencl.kernel_file_hierarchy_mmult_right", bpo::value(&opencl_kernel_file_hierarchy_mmult_right_)->default_value(opencl_kernel_file_hierarchy_mmult_right_), "Optional path to the hierarchy_mmult_right OpenCL kernel source code file.")
		("opencl.kernel_name_hierarchy_mmult_right", bpo::value(&opencl_kernel_name_hierarchy_mmult_right_)->default_value(opencl_kernel_name_hierarchy_mmult_right_), "Function name of the hierarchy_mmult_right OpenCL kernel.")
	;

	if (parse_config)
		parse(config_file_name);
}

} // namespace heom
