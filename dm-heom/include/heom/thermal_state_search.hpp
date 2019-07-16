// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_thermal_state_search_hpp
#define heom_thermal_state_search_hpp

#include "heom/instance.hpp"
#include <noma/ocl/helper.hpp>
#include <noma/ocl/kernel_wrapper.hpp>

namespace heom {

class thermal_state_search : public ocl::kernel_wrapper
{
public:
	thermal_state_search(ocl::helper& ocl, const boost::filesystem::path& thermal_state_search_file_name, const std::string& thermal_state_search_kernel_name, const std::string& source_header, const std::string& ocl_compile_options, const ocl::nd_range& range);
	thermal_state_search(ocl::helper& ocl, const std::string& source_header, const std::string& ocl_compile_options, const ocl::nd_range& range);

	void set_args(const instance& heom_instance, const complex_t epsilon, const complex_matrix_t& c, const complex_matrix_t& c_transposed, const complex_vector_t& eigenvalues);

	void run(cl::Buffer& d_mem_src, cl::Buffer& d_mem_dst);

private:
	static const std::string embedded_ocl_source_;
	static const std::string embedded_ocl_kernel_name_;

	cl::Buffer d_mem_c_;
	cl::Buffer d_mem_c_transposed_;
	cl::Buffer d_mem_eigenvalues_;
	cl::Buffer d_mem_pf_same_sum_;
};

} // namespace heom

#endif // heom_thermal_state_search_hpp