// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_communicator_hpp
#define heom_communicator_hpp

#include <numeric>

#include <mpi.h>
#include <noma/ocl/helper.hpp>
#include <noma/typa/matrix.hpp>

#include "debug.hpp"
#include "heom/hierarchy_partition.hpp"
#include "heom/instance.hpp"

namespace heom {

class instance;
class hierarchy_partition;

class communicator
{
public:
	communicator(const hierarchy_partition& partition, const instance& inst);
	~communicator();

	void neighbour_exchange(void* send_buffer, void* recv_buffer);

	size_t compute_size_byte() const { return compute_size_byte_; }
	size_t halo_size_byte() const { return halo_size_byte_; }

	void* send_buffer() const { return send_buffer_; }
	void* recv_buffer() const { return recv_buffer_; }

	const std::vector<rank_t>& sources() const { return sources_; }
	const std::vector<rank_t>& destinations() const { return destinations_; }

	const std::vector<size_t>& recv_mpi_type_sizes() const { return recv_mpi_type_sizes_; };
	const std::vector<size_t>& send_mpi_type_sizes() const { return send_mpi_type_sizes_; };

	void write_neighbour_exchange_table(std::ostream& os);
	void write_neighbour_exchange_table_raw(std::ostream& os);

	template<typename T> // bit-wise transferable type
	std::vector<T> gather(rank_t master_rank, const T& send_data);

	void barrier();

	static void init(int* argc, char*** argv);
	static void deinit();
	static int rank();
	static int ranks();

private:
	static bool initialised_;
	static int rank_; // rank of this instance's process
	static int ranks_; // number of ranks

	MPI_Comm neighbour_comm_; // neighbourhood communicator for this process

	std::vector<rank_t> sources_; // ranks to receive from into contiguous halo region of instance
	std::vector<rank_t> destinations_; // ranks to send to from non-contiguous compute region of instance

	std::vector<int> recv_counts_;
	std::vector<int> send_counts_;
	std::vector<MPI_Aint> recv_displacements_;
	std::vector<MPI_Aint> send_displacements_;
	std::vector<MPI_Datatype> recv_mpi_types_;
	std::vector<MPI_Datatype> send_mpi_types_;
	std::vector<size_t> recv_mpi_type_sizes_;
	std::vector<size_t> send_mpi_type_sizes_;

	size_t compute_size_byte_;
	size_t halo_size_byte_;

	void* send_buffer_; // start of compute region in instance's hierarchy buffer
	void* recv_buffer_; // start of halo region in instance's hierarchy buffer

	// MPI data types
	// NOTE: they are not actual types but data members from a language perspective
#ifdef HEOM_SINGLE_PRECISION
	const MPI_Datatype real_mpi_type_ = MPI_FLOAT;
#else
	const MPI_Datatype real_mpi_type_ = MPI_DOUBLE;
#endif

	MPI_Datatype hierarchy_node_mpi_type_; // one ado/sigma matrix

	parser::matrix<size_t> send_size_matrix_; // global send sizes, row i has all send sizes from rank i, column j the value for receiving rank j

	// error checker for mpi routines
	static void check_error(const int err, const std::string& msg);
};

template<typename T> // bit-wise transferable type
std::vector<T> communicator::gather(rank_t master_rank, const T& send_data)
{
	std::vector<T> recv_vec(ranks_);
	const size_t data_size = sizeof(T);
	check_error(MPI_Gather(reinterpret_cast<const void*>(&send_data),
	                       data_size,
	                       MPI_BYTE,
	                       reinterpret_cast<void*>(recv_vec.data()),
	                       data_size,
	                       MPI_BYTE,
	                       master_rank,
	                       MPI_COMM_WORLD),
	            "MPI_Gather");
	return recv_vec;
}

/**
 * Custom exception class for configuration errors.
 */
class communication_error : public std::runtime_error
{
public:
	communication_error(const std::string& msg) : runtime_error("heom::communication_error: " + msg) { };
};


/**
 * ODE pre-evaluate option adapter class for data exchange
 */
template<typename ODE_T>
class neighbour_exchange_action
{
public:
	neighbour_exchange_action(communicator& comm, instance& inst)
		: comm_(comm), inst_(inst)
	{
		DEBUG_ONLY( std::cout << "neighbour_exchange_action::neighbour_exchange_action(..): rank " << comm_.rank() << ": ctor called." << std::endl; )
	}

	void operator()(ODE_T& ode_inst, real_t time, real_t step_size, cl::Buffer& in, cl::Buffer& out)
	{
		DEBUG_ONLY( std::cout << "neighbour_exchange_action::operator()(..): rank " << comm_.rank() << ": called." << std::endl; )
		// NOTE: this is called pre-evaluate, i.e. in contains the last result
		assert(inst_.size_hierarchy_byte() == comm_.compute_size_byte() + comm_.halo_size_byte());

		void* send_buffer = comm_.send_buffer();
		void* recv_buffer = comm_.recv_buffer();

		// read hierarchy's compute nodes from in (device memory) into the instance's hierarchy buffer (host memory)
		// the halo nodes are written anyways during neighbour exchange
		{ // scope for timing
			bmt::timer timer;
			if (ode_inst.ocl_helper().zero_copy_available()) {
				send_buffer = ode_inst.ocl_helper().queue().enqueueMapBuffer(in, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, inst_.size_hierarchy_byte(), NULL, NULL, NULL);
				recv_buffer = reinterpret_cast<void*>(reinterpret_cast<std::int8_t*>(send_buffer) + comm_.compute_size_byte());
				// neighbour_exchange works directly in mapped memory now
			} else {
				ode_inst.ocl_helper().queue().enqueueReadBuffer(in, CL_TRUE, 0, comm_.compute_size_byte(), send_buffer, NULL, NULL);
			}
			buffer_read_stats_.add(timer);
		}

		{
			bmt::timer timer;
			comm_.neighbour_exchange(send_buffer, recv_buffer);
			neighbour_exchange_stats_.add(timer);
		}

		// write hierarchy's halo nodes back (with compute part size as offset)
		// compute nodes are unchanged after neighbour exchange
		// NOTE: offset (3. arg) is only for cl::Buffer and in byte
		{ // scope for timing
			bmt::timer timer;
			if (ode_inst.ocl_helper().zero_copy_available()) {
				ode_inst.ocl_helper().queue().enqueueUnmapMemObject(in, send_buffer, NULL, NULL);
				ode_inst.ocl_helper().queue().finish(); // make sure Unmap is through
			} else {
				ode_inst.ocl_helper().queue().enqueueWriteBuffer(in, CL_TRUE, comm_.compute_size_byte(), comm_.halo_size_byte(), recv_buffer, NULL, NULL);
			}
			buffer_write_stats_.add(timer);
		}
	}

	void write_runtime_stats(std::ostream& os) const
	{
		// header
		os << "name" << default_delimiter
		   << neighbour_exchange_stats_.header_string() << default_delimiter
		   << "rank" << default_delimiter
		   << "ranks" << default_delimiter
		   << "destinations" << default_delimiter
		   << "sources" << default_delimiter
		   << "send_size" << default_delimiter
		   << "recv_size" << '\n';

		std::stringstream additional_values;
		size_t send_size = std::accumulate(comm_.send_mpi_type_sizes().begin(), comm_.send_mpi_type_sizes().end(), 0);
		size_t recv_size = std::accumulate(comm_.recv_mpi_type_sizes().begin(), comm_.recv_mpi_type_sizes().end(), 0);;
		additional_values << comm_.rank() << default_delimiter
		                  << comm_.ranks() << default_delimiter
		                  << comm_.destinations().size() << default_delimiter
		                  << comm_.sources().size() << default_delimiter
		                  << send_size << default_delimiter
		                  << recv_size;

		// read/map
		os << "buffer_read_stats" << default_delimiter
		   << buffer_read_stats_.string() << default_delimiter
		   << additional_values.str() << '\n';

		// kernel_heom_ode_stats
		os << "neighour_exchange_stats" << default_delimiter
		   << neighbour_exchange_stats_.string() << default_delimiter
		   << additional_values.str() << '\n';

		// write/unmap stats
		os << "buffer_write_stats" << default_delimiter
		   << buffer_write_stats_.string() << default_delimiter
		   << additional_values.str() << std::endl;
	}

	void write_runtime_summary(std::ostream& os)
	{
		write_time_statistics_summary(buffer_read_stats_, "read buffer", os);
		write_time_statistics_summary(neighbour_exchange_stats_, "neighbour exchange", os);
		write_time_statistics_summary(buffer_write_stats_, "write buffer", os);
	}

	statistics_vector stats() const
	{
		statistics_vector stats_vec;

		stats_vec.push_back(statistics_vector::value_type(buffer_read_stats_, "read buffer"));
		stats_vec.push_back(statistics_vector::value_type(neighbour_exchange_stats_, "neighbour exchange"));
		stats_vec.push_back(statistics_vector::value_type(buffer_write_stats_, "write buffer"));

		return stats_vec;
	}

private:
	communicator& comm_;
	instance& inst_;
	bmt::statistics buffer_read_stats_;
	bmt::statistics neighbour_exchange_stats_;
	bmt::statistics buffer_write_stats_;
};

} // namespace heom

#endif // heom_communicator_hpp
