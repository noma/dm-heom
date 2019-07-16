// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_make_observer_hpp
#define heom_make_observer_hpp

#include <memory>

#include "heom/matrix_complete_observer.hpp"
#include "heom/matrix_diagonal_observer.hpp"
#include "heom/matrix_trace_observer.hpp"
#include "heom/observer.hpp"
#include "heom/observation_type.hpp"
#include "heom/hierarchy_norm_observer.hpp"
#include "heom/hierarchy_complete_observer.hpp"

namespace heom {

/**
 * Create C++ observer type according to observation_type value
 * @tparam SOLVER_T solver type to observe, must implement solver concept
 * @tparam Args constructor argument types (NOTE: all created types must still have the same constructor signature)
 * @param observation_type determines the C++ observer type to create
 * @param args constructor arguments
 * @return std::unique_ptr<observer> to the newly created observer instance
 */
template<typename SOLVER_T, typename... Args>
std::unique_ptr<observer> make_unique_observer(const observation_type_t& observation_type, Args&&... args)
{

	switch (observation_type) {
	case observation_type_t::matrix_complete:
		return std::unique_ptr<observer> { new matrix_complete_observer<SOLVER_T>(std::forward<Args>(args)...) };
	case observation_type_t::matrix_diagonal:
		return std::unique_ptr<observer> { new matrix_diagonal_observer<SOLVER_T>(std::forward<Args>(args)...) };
	case observation_type_t::matrix_trace_id:
		return std::unique_ptr<observer> { new matrix_trace_observer<SOLVER_T>(std::forward<Args>(args)...) };
	case observation_type_t::hierarchy_norm:
		return std::unique_ptr<observer> { new hierarchy_norm_observer<SOLVER_T>(std::forward<Args>(args)...) };
	case observation_type_t::hierarchy_complete:
		return std::unique_ptr<observer> { new hierarchy_complete_observer<SOLVER_T>(std::forward<Args>(args)...) };
	case observation_type_t::matrix_trace_transient_absorption:
		throw std::runtime_error("make_unique_observer(): error: matrix_trace_transient_absorption cannot be created using make_unique_observer.");
	default:
		throw std::runtime_error("make_unique_observer(): error: encountered unknown observation type while setting up observers.");
	}
}

} // namespace heom

#endif // heom_make_observer_hpp
