// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_make_file_observer_list_hpp
#define heom_make_file_observer_list_hpp

#include <memory>
#include <vector>

#include "heom/config.hpp"
#include "heom/file_observer_wrapper.hpp"
#include "heom/instance.hpp"
#include "heom/make_observer.hpp"
#include "heom/observer_list.hpp"

namespace heom {

template<typename SOLVER_T>
observer_list make_file_observer_list(const config::obervation_type_list& observations, const hierarchy_graph& graph, SOLVER_T& solver)
{
	observer_list obs_list(observations.get().size());
	for (auto& pair : observations.get()) {
		auto& obs_type = pair.get().first;
		auto& filename = pair.get().second;

		// specific types, need specific constructor arguments
		switch (obs_type) {
			// NOTE: all the same ctor taking a matrix dimension, i.e. system_sites()
			case observation_type_t::matrix_complete:
			case observation_type_t::matrix_diagonal:
			case observation_type_t::matrix_trace_id:
			case observation_type_t::hierarchy_norm:
			case observation_type_t::hierarchy_complete:
				obs_list.add(std::unique_ptr<observer> { new file_observer_wrapper(make_unique_observer<SOLVER_T>(obs_type, graph, solver), filename) });
				break;
			case observation_type_t::matrix_trace_transient_absorption:
				throw std::runtime_error("make_unique_observer(): error: matrix_trace_transient_absorption cannot be created using make_file_observer_list.");
			default:
				throw std::runtime_error("make_observer_list(): error: encountered unknown observation type while setting up observers.");
		}
	}

	return obs_list;
}

template<typename SOLVER_T>
observer_list make_file_observer_list(const config& cfg, const hierarchy_graph& graph, SOLVER_T& solver)
{
	return make_file_observer_list(cfg.observations(), graph, solver);
}


} // namespace heom

#endif // heom_make_file_observer_list_hpp
