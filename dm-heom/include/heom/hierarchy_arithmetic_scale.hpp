// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
// Copyright (c) 2016-2017 Lucas Deecke, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_hierarchy_arithmetic_scale_hpp
#define heom_hierarchy_arithmetic_scale_hpp

#include <noma/ocl/helper.hpp>
#include <noma/ocl/kernel_wrapper.hpp>

#include "heom/instance.hpp"

namespace heom {

class hierarchy_arithmetic_scale : public ocl::kernel_wrapper
{
public:
	hierarchy_arithmetic_scale(ocl::helper& ocl, const boost::filesystem::path& hierarchy_arithmetic_scale_file_name, const std::string& hierarchy_arithmetic_scale_kernel_name, const std::string& source_header, const std::string& ocl_compile_options, const ocl::nd_range& range);
	hierarchy_arithmetic_scale(ocl::helper& ocl, const std::string& source_header, const std::string& ocl_compile_options, const ocl::nd_range& range);

	void run(cl::Buffer& d_mem_a, complex_t b, cl::Buffer& d_mem_dst);

private:
	static const std::string embedded_ocl_source_;
	static const std::string embedded_ocl_kernel_name_;
};

} // namespace heom

#endif // heom_hierarchy_arithmetic_scale_hpp
