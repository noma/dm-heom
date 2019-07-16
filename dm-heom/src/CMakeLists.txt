# This file is part of DM-HEOM (https://github.com/noma/dm-heom)
#
# Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
#
# Licensed under the 3-clause BSD License, see accompanying LICENSE,
# CONTRIBUTORS.md, and README.md for further information.

message(STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR})
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(heom)

### Non-distributed Applications
# population dynamics
add_executable(app_population_dynamics app_population_dynamics.cpp)
target_link_libraries(app_population_dynamics dm-heom)

# rate kernel
add_executable(app_rate_kernel app_rate_kernel.cpp)
target_link_libraries(app_rate_kernel dm-heom)

# linear absorption
add_executable(app_linear_absorption app_linear_absorption.cpp)
target_link_libraries(app_linear_absorption dm-heom)

# circular dichroism
add_executable(app_circular_dichroism app_circular_dichroism.cpp)
target_link_libraries(app_circular_dichroism dm-heom)

# transient absorption
add_executable(app_transient_absorption app_transient_absorption.cpp)
target_link_libraries(app_transient_absorption dm-heom)

# thermal state search
add_executable(app_thermal_state_search app_thermal_state_search.cpp)
target_link_libraries(app_thermal_state_search dm-heom)

# static fluorescence
add_executable(app_static_fluorescence app_static_fluorescence.cpp)
target_link_libraries(app_static_fluorescence dm-heom)

# two dimensional spectra
add_executable(app_two_dimensional_spectra app_two_dimensional_spectra.cpp)
target_link_libraries(app_two_dimensional_spectra dm-heom)

### Utilities/Tools

add_executable(util_config_help util_config_help.cpp)
target_link_libraries(util_config_help dm-heom)

add_executable(util_config_check util_config_check.cpp)
target_link_libraries(util_config_check dm-heom)

add_executable(util_graph_to_file util_graph_to_file.cpp)
target_link_libraries(util_graph_to_file dm-heom)

add_executable(util_partitioning_info util_partitioning_info.cpp)
target_link_libraries(util_partitioning_info dm-heom)

### Distributed applications

if (HEOM_ENABLE_MPI)
	# population dynamics
	add_executable(app_mpi_population_dynamics app_mpi_population_dynamics.cpp)
	target_link_libraries(app_mpi_population_dynamics dm-heom_mpi)
	if (MPI_CXX_COMPILE_FLAGS)
		set_target_properties(app_mpi_population_dynamics PROPERTIES COMPILE_FLAGS ${MPI_CXX_COMPILE_FLAGS})
	endif ()
	if (MPI_CXX_LINK_FLAGS)
		set_target_properties(app_mpi_population_dynamics PROPERTIES LINK_FLAGS ${MPI_CXX_LINK_FLAGS})
	endif ()
endif ()


### Tests

if (NOT DEFINED HEOM_BUILD_TESTS)
	if (NOT(CMAKE_BUILD_TYPE STREQUAL "Release")) # disable in Release mode by default
		set(HEOM_BUILD_TESTS TRUE)
	else ()
		set(HEOM_BUILD_TESTS FALSE)
	endif ()
endif ()

if (HEOM_BUILD_TESTS)
	message(STATUS "Building with tests")

	if (HEOM_ENABLE_MPI)
	endif ()

endif()
