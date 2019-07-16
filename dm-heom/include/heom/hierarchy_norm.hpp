// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_hierarchy_norm_hpp
#define heom_hierarchy_norm_hpp

#include <noma/ocl/helper.hpp>
#include <noma/ocl/kernel_wrapper.hpp>

#include "heom/common.hpp"

namespace heom {

class hierarchy_norm : public ocl::kernel_wrapper
{
public:
	hierarchy_norm(ocl::helper& ocl, const std::string& source_header, const std::string& ocl_compile_options, const ocl::nd_range& range);
	hierarchy_norm(ocl::helper& ocl, const boost::filesystem::path& hierarchy_norm_file_name, const std::string& hierarchy_norm_kernel_name, const std::string& source_header, const std::string& ocl_compile_options, const ocl::nd_range& range);

	void update(real_t threshold, cl::Buffer& d_mem_norm, cl::Buffer& d_mem_mask, cl::Buffer& hierarchy);
	void update(cl::Buffer& d_mem_norm, cl::Buffer& hierarchy);

private:
	static const std::string embedded_ocl_source_;
	static const std::string embedded_ocl_kernel_name_;
};

} // namespace heom

#endif // heom_hierarchy_norm_hpp
