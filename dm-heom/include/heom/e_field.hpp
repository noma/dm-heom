// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_e_field_hpp
#define heom_e_field_hpp

#include <valarray>

#include "heom/common.hpp"

// forward declaration instead of header
namespace cl {
	class Buffer;
}

namespace heom {

/**
 * Representation of a time-dependent electric field, used to modify the
 * Hamiltonian in pump/probe laser pulse scenarios.
 *
 * Implemented as a functor (callable) that meets the concept of an
 * pre_evaluation_action of the ODE. It is called with the exact time
 * (time + step_size) with which the ODE is going to be evaluated just before
 * its evaluation to update the Hamiltonian with the e-field at that same
 * point in time.
 */
template<typename CONFIG_T, typename ODE_T>
class e_field {
public:
	e_field(CONFIG_T& config, size_t phase_index, const complex_matrix_t& dipole_plus, const complex_matrix_t& dipole_minus, const complex_matrix_t& hamiltonian)
	  : config_(config),
	    size_(dipole_plus.rows() * dipole_plus.cols()),
	    phase_index_(phase_index),
	    d_plus_(dipole_plus.data(), size_),
	    d_minus_(dipole_minus.data(), size_),
	    hamiltonian_in_(hamiltonian.data(), size_)
	{}

	// compute time dependent e-field
	complex_t compute_field(real_t t) const
	{
		// compute field from pump/probe
		// EEpump[t_] := Epump0*Exp[-(t-tpulse)*(t-tpulse)/(2*tpump*tpump)];
		// EEprobe[t_, tdel_] := Eprobe0*Exp[-(t-tpulse-tdel)*(t-tpulse-tdel)/(2*tprobe*tprobe)];
		complex_t field =
			config_.probe_pulse_efield() * exp(-pow((t - config_.probe_pulse_time_center()), 2)
			                                   / (2.0 * pow(config_.probe_pulse_time_full_width(), 2)))
			* exp(complex_t(0, 1) * (config_.probe_pulse_phases()[phase_index_] + config_.probe_pulse_frequency() * t))
			+
			config_.pump_pulse_efield() * exp(-pow((t - config_.pump_pulse_time_center()), 2)
			                                  / (2.0 * pow(config_.pump_pulse_time_full_width(), 2)))
			* exp(complex_t(0, 1) * (config_.pump_pulse_phases()[phase_index_] + config_.pump_pulse_frequency() * t));

		return field;
	}

	// ode pre-evaluate-action compatible interface
	void operator()(ODE_T& ode, real_t time, real_t step_size, cl::Buffer& in, cl::Buffer& out) const
	{
		complex_t field = compute_field(time + step_size);
		// compute new hamiltonian
		// TODO: gcc 5.2.1: using auto as result type here segfaults
		// segfault:
//		auto hamiltonian_out = hamiltonian_in + field * d_plus + std::conj(field) * d_minus; // segfault
		// segfault:
//		auto x = field * d_plus + std::conj(field) * d_minus;
//		auto hamiltonian_out = hamiltonian_in + x;
		// works:
//		auto x = field * d_plus;
//		auto y = std::conj(field) * d_minus;
//		auto hamiltonian_out = hamiltonian_in + x;
		// works:
		std::valarray<complex_t> hamiltonian_out = hamiltonian_in_ - field * d_minus_ - std::conj(field) * d_plus_;

		// update hamiltonian in instance
		ode.get_heom_instance().set_hamiltonian(hamiltonian_out);

		// notify ode about change
		ode.update_hamiltonian();
	}

private:
	const CONFIG_T& config_;
	const size_t size_;
	const size_t phase_index_;
	const std::valarray<complex_t> d_plus_;
	const std::valarray<complex_t> d_minus_;
	const std::valarray<complex_t> hamiltonian_in_;
};

} // namespace heom

#endif // heom_e_field_hpp

