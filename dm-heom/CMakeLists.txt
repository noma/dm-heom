# This file is part of DM-HEOM (https://github.com/noma/dm-heom)
#
# Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
#
# Licensed under the 3-clause BSD License, see accompanying LICENSE,
# CONTRIBUTORS.md, and README.md for further information.

cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR) # verfified 2017-08-09

project(HEOM)

# Set default build type to "Release"
set(HEOM_DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "Setting build type to '${HEOM_DEFAULT_BUILD_TYPE}' as none was specified.")
	set(CMAKE_BUILD_TYPE "${HEOM_DEFAULT_BUILD_TYPE}" CACHE
			STRING "Choose the type of build." FORCE)
	# Set the possible values of build type for cmake-gui
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
			"Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# test compiler versions
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9.0")
		message(FATAL_ERROR "GCC version 4.9.0 or higher required.")
	endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.3.0")
		message(FATAL_ERROR "Clang version 3.3.0 or higher required.")
	endif()
	execute_process( 
		COMMAND bash -c "echo -n $(dirname $(which g++))/.."
		OUTPUT_VARIABLE GCC_PATH
	)
	message(WARNING "Clang detected, setting '--gcc-toolchain=${GCC_PATH}' option, make sure version is >=4.9.0.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --gcc-toolchain=${GCC_PATH}")
else()
	message(WARNING "Using untested compiler, full C++11 compliance including the standard library is required. When using LLVM/Clang or Intel with a GNU standard library, make sure the used GCC installation is at least version 4.9.0")
endif()

# location of all thirdparty components
set(THIRDPARTY_DIR ${PROJECT_SOURCE_DIR}/../thirdparty)

# location of OpenCL kernels
set(OpenCL_KERNEL_DIR ${PROJECT_SOURCE_DIR}/../cl)

# search for additional cmake modules
# git version hash
list(APPEND CMAKE_MODULE_PATH "${THIRDPARTY_DIR}/cmake/rpavlik/")

# set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
message(STATUS "CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR})

# compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
if (NOT(CMAKE_BUILD_TYPE STREQUAL "Release")) # enable debug symbols for all but Release build
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

# Boost
#set(Boost_USE_STATIC_LIBS OFF)
#set(Boost_USE_MULTITHREADED ON)
#set(Boost_USE_STATIC_RUNTIME OFF)

if (USE_BOOST_ROOT)
	message(FATAL_ERROR "USE_BOOST_ROOT is deprecated and now and no longer supported.")
endif ()

find_package(Boost REQUIRED COMPONENTS program_options filesystem)
# avoid compiler warnings from Boost.System for unused symbols
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_SYSTEM_NO_DEPRECATED=1")

message(STATUS "Boost_INCLUDE_DIRS: " ${Boost_INCLUDE_DIRS})
message(STATUS "Boost_LIBRARY_DIRS: " ${Boost_LIBRARY_DIRS})
message(STATUS "Boost_LIBRARIES: " ${Boost_LIBRARIES})

if (HEOM_ENABLE_MPI)
	find_package(MPI REQUIRED)
	include_directories(${MPI_CXX_INCLUDE_PATH})
endif ()

# Eigen
set(Eigen_INCLUDE_DIRS ../thirdparty/eigen/)
add_library(eigen INTERFACE)
target_include_directories (eigen INTERFACE ${Eigen_INCLUDE_DIRS})

# project wide include dirs
include_directories(${PROJECT_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS})

# thirdparty CMake projects
# misc. helpers
add_subdirectory(../thirdparty/misc/ ${CMAKE_CURRENT_BINARY_DIR}/build.noma_misc)
# benchmark timer
add_subdirectory(../thirdparty/bmt/ ${CMAKE_CURRENT_BINARY_DIR}/build.noma_bmt)
# type parser library
add_subdirectory(../thirdparty/typa/ ${CMAKE_CURRENT_BINARY_DIR}/build.noma_typa)
# OpenCL helper library
add_subdirectory(../thirdparty/ocl/ ${CMAKE_CURRENT_BINARY_DIR}/build.noma_ocl)
# numerics library
add_subdirectory(../thirdparty/num/ ${CMAKE_CURRENT_BINARY_DIR}/build.noma_num)

# see CMakeLists.txt in src
add_subdirectory(src)

message(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})
