# This file is part of DM-HEOM (https://github.com/noma/dm-heom)
#
# Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
#
# Licensed under the 3-clause BSD License, see accompanying LICENSE,
# CONTRIBUTORS.md, and README.md for further information.

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(libdm-heom LANGUAGES CXX)

message(STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR})
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR})

# git version source file (avoids recompilation as compared to a define)
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
#set(GIT_VERSION "${GIT_REFSPEC}")
set(GIT_VERSION "${GIT_SHA1}")
git_describe(GIT_DESCRIBE)
#set(GIT_VERSION "${GIT_DESCRIBE}")
git_get_exact_tag(GIT_EXACT_TAG)
#set(GIT_VERSION "${GIT_EXACT_TAG}")
git_local_changes(GIT_LOCAL_CHANGES)
set(GIT_LOCAL_CHANGES "${GIT_LOCAL_CHANGES}")
# generate git version source file from template
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/git_version.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/git_version.cpp" @ONLY)

# generate compiler version header (usese CMAKE_CXX_COMPILER_ID and CMAKE_CXX_COMPILER_VERSION)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/compiler_version.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/compiler_version.cpp" @ONLY)

# all library source files in alphabetical order
set(HEOM_LIB_SRC
		${CMAKE_CURRENT_BINARY_DIR}/compiler_version.cpp ${CMAKE_CURRENT_BINARY_DIR}/git_version.cpp
		circular_dichroism_config.cpp
		command_line.cpp
		common.cpp
		config.cpp
		config_base.cpp
		constants.cpp
		dipole_matrix.cpp
		eigen_decomposition.cpp
		error_norm.cpp
		filtering_strategy.cpp
		handle_main_exception.cpp
		hierarchy_arithmetic_add.cpp
		hierarchy_arithmetic_dot_product.cpp
		hierarchy_arithmetic_mmult_left.cpp
		hierarchy_arithmetic_mmult_right.cpp
		hierarchy_arithmetic_scale.cpp
		hierarchy_arithmetic_subtract.cpp
		hierarchy_graph.cpp
		hierarchy_mask.cpp
		hierarchy_norm.cpp
		hierarchy_partition.cpp
		instance.cpp
		integer_partition.cpp
		observation_type.cpp
		observer.cpp
		ocl_config.cpp
		ode.cpp
		population_dynamics_solver.cpp
		pulse_shape.cpp
		program_task.cpp
		rk_weighted_add.cpp
		sites_to_states.cpp
		spectra_pathway.cpp
		thermal_state_search.cpp
		thermal_state_search_solver.cpp

		dipole_config_entries.cpp
		dipole_matrix_config_entries.cpp
		laser_config_entries.cpp
		linear_absorption_config.cpp
		population_dynamics_config.cpp
		rate_kernel_config.cpp
		static_fluorescence_config.cpp
		thermal_state_search_config.cpp
		transient_absorption_config.cpp
		two_dimensional_spectra_config.cpp
		)

set(HEOM_LIB_MPI_SRC
		communicator.cpp
		distributed_command_line.cpp
		)


# library depencies
set(HEOM_LIB_DEPS noma_debug noma_memory noma_bmt noma_typa noma_ocl noma_num eigen)

# OpenCL kernels are compiled into the application
set(OpenCL_KERNEL_DIR ${PROJECT_SOURCE_DIR}/../../cl/heom) # where the kernel files are, convenience
set(OpenCL_KERNEL_HEADER_DIR ${CMAKE_CURRENT_BINARY_DIR}/cl) # where the generated header file will be put, used by create_opencl_kernel_header
file(MAKE_DIRECTORY ${OpenCL_KERNEL_HEADER_DIR})
include_directories(${OpenCL_KERNEL_HEADER_DIR})

# heom_ode.cl used by ode.hpp/cpp
create_opencl_kernel_header(${OpenCL_KERNEL_DIR}/heom_ode.cl ${OpenCL_KERNEL_HEADER_DIR} KERNEL_HEADER_heom_ode)

# thermal_state_search.cl used by thermal_state_search.hpp/cpp
create_opencl_kernel_header(${OpenCL_KERNEL_DIR}/thermal_state_search.cl ${OpenCL_KERNEL_HEADER_DIR} KERNEL_HEADER_thermal_state_search)

# hierarchy_norm.cl used by hierarchy_norm.hpp/cpp
create_opencl_kernel_header(${OpenCL_KERNEL_DIR}/hierarchy_norm.cl ${OpenCL_KERNEL_HEADER_DIR} KERNEL_HEADER_hierarchy_norm)

# error_norm.cl used by error_norm.hpp/cpp
create_opencl_kernel_header(${OpenCL_KERNEL_DIR}/error_norm.cl ${OpenCL_KERNEL_HEADER_DIR} KERNEL_HEADER_error_norm)

# hierarchy_arithmetic_add.cl used by hierarchy_arithmetic_add.hpp/cpp
create_opencl_kernel_header(${OpenCL_KERNEL_DIR}/hierarchy_arithmetic_add.cl ${OpenCL_KERNEL_HEADER_DIR} KERNEL_HEADER_hierarchy_arithmetic_add)

# hierarchy_arithmetic_subtract.cl used by hierarchy_arithmetic_subtract.hpp/cpp
create_opencl_kernel_header(${OpenCL_KERNEL_DIR}/hierarchy_arithmetic_subtract.cl ${OpenCL_KERNEL_HEADER_DIR} KERNEL_HEADER_hierarchy_arithmetic_subtract)

# hierarchy_arithmetic_scale.cl used by hierarchy_arithmetic_scale.hpp/cpp
create_opencl_kernel_header(${OpenCL_KERNEL_DIR}/hierarchy_arithmetic_scale.cl ${OpenCL_KERNEL_HEADER_DIR} KERNEL_HEADER_hierarchy_arithmetic_scale)

# hierarchy_arithmetic_dot_product.cl used by hierarchy_arithmetic_dot_product.hpp/cpp
create_opencl_kernel_header(${OpenCL_KERNEL_DIR}/hierarchy_arithmetic_dot_product.cl ${OpenCL_KERNEL_HEADER_DIR} KERNEL_HEADER_hierarchy_arithmetic_dot_product)

# hierarchy_arithmetic_mmult_left.cl used by hierarchy_arithmetic_mmult_left.hpp/cpp
create_opencl_kernel_header(${OpenCL_KERNEL_DIR}/hierarchy_arithmetic_mmult_left.cl ${OpenCL_KERNEL_HEADER_DIR} KERNEL_HEADER_hierarchy_arithmetic_mmult_left)

# hierarchy_arithmetic_mmult_right.cl used by hierarchy_arithmetic_mmult_right.hpp/cpp
create_opencl_kernel_header(${OpenCL_KERNEL_DIR}/hierarchy_arithmetic_mmult_right.cl ${OpenCL_KERNEL_HEADER_DIR} KERNEL_HEADER_hierarchy_arithmetic_mmult_right)

# rk_weighted_add.cl used by hierarchy_arithmetic_mmult_right.hpp/cpp
create_opencl_kernel_header(${OpenCL_KERNEL_DIR}/rk_weighted_add.cl ${OpenCL_KERNEL_HEADER_DIR} KERNEL_HEADER_rk_weighted_add)

# NOTE: one meta target for OpenCL_Kernel_Headers the applications to depend on avoid conflicts due to multiple generation attempts in parallel build (for every object for every app). The dependencies between source file and included generated header are figured out by CMake, so set_property with OBJECT_DEPENDENCY is needed.
add_custom_target(OpenCL_Kernel_Headers DEPENDS ${OpenCL_KERNEL_HEADER_DIR} ${KERNEL_HEADER_heom_ode} ${KERNEL_HEADER_thermal_state_search} ${KERNEL_HEADER_hierarchy_norm} ${KERNEL_HEADER_error_norm} ${KERNEL_HEADER_hierarchy_arithmetic_add} ${KERNEL_HEADER_hierarchy_arithmetic_subtract} ${KERNEL_HEADER_hierarchy_arithmetic_scale} ${KERNEL_HEADER_hierarchy_arithmetic_dot_product} ${KERNEL_HEADER_hierarchy_arithmetic_mmult_left} ${KERNEL_HEADER_hierarchy_arithmetic_mmult_right} ${KERNEL_HEADER_rk_weighted_add})

### HEOM library

# static dm-heom library (non-distributed)
add_library(dm-heom STATIC ${HEOM_LIB_SRC})
target_link_libraries(dm-heom PUBLIC ${Boost_LIBRARIES} ${HEOM_LIB_DEPS})
add_dependencies(dm-heom OpenCL_Kernel_Headers)
target_include_directories(dm-heom PUBLIC include)

set_target_properties(dm-heom PROPERTIES
		CXX_STANDARD 11
		CXX_STANDARD_REQUIRED YES
		CXX_EXTENSIONS NO
		)

# static dm-heom MPI commmunication library
if (HEOM_ENABLE_MPI)
	add_library(dm-heom_mpi STATIC ${HEOM_LIB_MPI_SRC})
	target_link_libraries(dm-heom_mpi PUBLIC dm-heom ${MPI_CXX_LIBRARIES})
	target_include_directories(dm-heom_mpi PUBLIC include)

	set_target_properties(dm-heom_mpi PROPERTIES
			CXX_STANDARD 11
			CXX_STANDARD_REQUIRED YES
			CXX_EXTENSIONS NO
			)
endif()
