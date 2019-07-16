// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/communicator.hpp"

#include "heom/instance.hpp"
#include "heom/hierarchy_partition.hpp"


#include "debug.hpp"

namespace heom {

const std::string mpi_error_to_string(int err)
{
	switch (err)
	{
		case MPI_SUCCESS:
			return "No error (MPI_SUCCESS).";
		case MPI_ERR_BUFFER:
			return "Invalid buffer pointer (MPI_ERR_BUFFER).";
		case MPI_ERR_COUNT:
			return "Invalid count argument (MPI_ERR_COUNT).";
		case MPI_ERR_TYPE:
			return "Invalid datatype argument (MPI_ERR_TYPE).";
		case MPI_ERR_TAG:
			return "Invalid tag argument (MPI_ERR_TAG).";
		case MPI_ERR_COMM:
			return "Invalid communicator (MPI_ERR_COMM).";
		case MPI_ERR_RANK:
			return "Invalid rank (MPI_ERR_RANK).";
		case MPI_ERR_REQUEST:
			return "Invalid request (handle) (MPI_ERR_REQUEST).";
		case MPI_ERR_ROOT:
			return "Invalid root (MPI_ERR_ROOT).";
		case MPI_ERR_GROUP:
			return "Invalid group (MPI_ERR_GROUP).";
		case MPI_ERR_OP:
			return "Invalid operation (MPI_ERR_OP).";
		case MPI_ERR_TOPOLOGY:
			return "Invalid topology (MPI_ERR_TOPOLOGY).";
		case MPI_ERR_DIMS:
			return "Invalid dimension argument (MPI_ERR_DIMS).";
		case MPI_ERR_ARG:
			return "Invalid argument of some other kind (MPI_ERR_ARG).";
		case MPI_ERR_UNKNOWN:
			return "Unknown error (MPI_ERR_UNKNOWN).";
		case MPI_ERR_TRUNCATE:
			return " (Message truncated on receive).";
		case MPI_ERR_OTHER:
			return "Known error not in this list (MPI_ERR_OTHER).";
		case MPI_ERR_INTERN:
			return "Internal MPI (implementation) error (MPI_ERR_INTERN).";
		case MPI_ERR_IN_STATUS:
			return "Error code is in status (MPI_ERR_IN_STATUS).";
		case MPI_ERR_PENDING:
			return "Pending request (MPI_ERR_PENDING).";
		case MPI_ERR_LASTCODE:
			return "Last error code (MPI_ERR_LASTCODE).";
		default:
			return "Unknown error code.";
	}
}

bool communicator::initialised_ = false;
int communicator::rank_ = -1;
int communicator::ranks_ = -1;

void communicator::check_error(const int err, const std::string& msg)
{
	if (err != MPI_SUCCESS)
		throw communication_error(msg + mpi_error_to_string(err));
}

void communicator::init(int* argc, char*** argv)
{
	check_error(MPI_Init(argc, argv), "MPI_Init");

	check_error(MPI_Comm_rank(MPI_COMM_WORLD, &rank_), "MPI_Comm_rank");
	check_error(MPI_Comm_size(MPI_COMM_WORLD, &ranks_), "MPI_Comm_size");

	initialised_ = true;
}

void communicator::deinit()
{
	if (!initialised_)
		throw communication_error("communicator::deinit() called without prior init().");

	//check_error(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier"); // NOTE: activate for debugging, with output before and after
	check_error(MPI_Finalize(), "MPI_Finalize");
}

int communicator::rank()
{
	return rank_;
}

int communicator::ranks()
{
	return ranks_;
}


// NOTE: the idea here is to
// - build send/recv data types so MPI handles packing/unpacking,
// - setup neighbourhood communicators to reflect the partition graph, and
// - then use MPI_Neighbour_alltoallw on those communictor with the data types to do the neighbour exchange
// NOTE: out = sent by this process, in = received by this process
communicator::communicator(const hierarchy_partition& partition, const instance& inst)
	: sources_(partition.in_neighbour_ranks().begin(), partition.in_neighbour_ranks().end()), // collect source ranks into a vector that can be used in MPI-calls via data()
	  destinations_(partition.out_neighbour_ranks().begin(), partition.out_neighbour_ranks().end()), // same for destination rank, see above
	  compute_size_byte_(partition.compute_nodes() * inst.size_hierarchy_top_byte()),
	  halo_size_byte_(partition.halo_nodes() * inst.size_hierarchy_top_byte()),
	  send_buffer_(reinterpret_cast<void*>(inst.hierarchy())),
	  recv_buffer_(reinterpret_cast<void*>(reinterpret_cast<std::int8_t*>(inst.hierarchy()) + compute_size_byte_)), // NOTE: cast to int8_t for byte-sized pointer arithmetic
	  send_size_matrix_(ranks_, ranks_, 0)
{
	// set container sizes to avoid re-allocations
	recv_counts_.reserve(sources_.size());
	send_counts_.reserve(destinations_.size());
	recv_displacements_.reserve(sources_.size());
	send_displacements_.reserve(destinations_.size());
	recv_mpi_types_.reserve(sources_.size());
	send_mpi_types_.reserve(destinations_.size());
	recv_mpi_type_sizes_.reserve(sources_.size());
	send_mpi_type_sizes_.reserve(destinations_.size());

	// set up neighbourhood communicator from partition graph
	check_error(MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD, // comm_old
	                                           sources_.size(), // indegree
	                                           sources_.data(), // sources[]
	                                           reinterpret_cast<const int *>(MPI_UNWEIGHTED), // sourceweights, NOTE: in some implementations, MPI_UNWEIGHTED has the wrong type
	                                           destinations_.size(), // outdegree
	                                           destinations_.data(), // destinations[]
	                                           reinterpret_cast<const int *>(MPI_UNWEIGHTED), // destweights, NOTE: in some implementations, MPI_UNWEIGHTED has the wrong type
	                                           MPI_INFO_NULL, // info
	                                           0, // reorder
	                                           &neighbour_comm_), // comm_dist_graph
	            "MPI_Dist_graph_create_adjacent");

	// create MPI type for a single hierarchy node, aka sigma matrix
	check_error(MPI_Type_contiguous(inst.size_hierarchy_top_byte() / sizeof(real_t), real_mpi_type_, &hierarchy_node_mpi_type_), "MPI_Type_contiguous");
	check_error(MPI_Type_commit(&hierarchy_node_mpi_type_), "MPI_Type_commit");

	DEBUG_ONLY(
		int type_size = 0;
		check_error(MPI_Type_size(real_mpi_type_, &type_size), "MPI_Type_size");
		assert(static_cast<size_t>(type_size) == sizeof(real_t));
		check_error(MPI_Type_size(hierarchy_node_mpi_type_, &type_size), "MPI_Type_size");
		assert(static_cast<size_t>(type_size) == inst.size_hierarchy_top_byte());
	)

	// setup send/out/destinations stuff (complicated, non-contiguous sub-sets of nodes computed on this rank, needed by others, matched by a contiguous receive data type on another rank)

	// create MPI data types for each rank this one is sending to
	ASSERT_ONLY( rank_t last_destination = -1; )
	for (const auto& destination : destinations_)
	{
		assert(destination > last_destination); // make sure the expected rank order is valid, i.e. numerical order
		// collect blocks of of contiguous hierarchy nodes for that each neighbour rank

		const auto& send_uids = partition.rank_to_out_neighbours().at(destination); // set of nodes this rank needs to send to destination, global uids

		int current_block_length = 0;
		std::vector<int> displacements; // in hierarchy nodes
		std::vector<int> block_lengths; // in hierarchy nodes

		uid_t last_lid = hierarchy_graph::invalid_uid; // track last uid to detect contiguously aligned blocks of local hierarchy nodes
		for (const auto& current_uid : send_uids)
		{
			rank_t current_lid = partition.uid_to_lid().at(current_uid);
			assert(current_lid < partition.compute_nodes()); // must be a locally computed hierarchy node
			if (current_lid != (last_lid + 1) || last_lid == hierarchy_graph::invalid_uid) // block start/end non-contiguous lid sequence, or first block
			{
				// end block, if there was one
				if (last_lid != hierarchy_graph::invalid_uid)
				{
					assert(current_block_length > 0);
					block_lengths.push_back(current_block_length); // add the block lengths
				}

				// start new block
				current_block_length = 0; // reset block length
				displacements.push_back(current_lid); // block starts at displacement
			}
			last_lid = current_lid;
			++current_block_length; // increment block length in any case
		}
		block_lengths.push_back(current_block_length); // finish last block
		block_lengths.shrink_to_fit();
		displacements.shrink_to_fit();
		assert(block_lengths.size() == displacements.size());

		send_counts_.push_back(1); // NOTE: currently exactly one instance per type is required
		send_displacements_.push_back(0); // NOTE: currently not needed, first block of MPI_Type_indexed contains a displacement already, i.e. each type is relative to the start of send_buffer

		MPI_Datatype new_mpi_type;
		check_error(MPI_Type_indexed(block_lengths.size(), // int count
		                             block_lengths.data(), // int* array_of_blocklengths,
		                             displacements.data(), // int *array_of_displacements,
		                             hierarchy_node_mpi_type_, // MPI_Datatype oldtype,
		                             &new_mpi_type), // MPI_Datatype *newtype)
		            "MPI_Type_indexed");
		check_error(MPI_Type_commit(&new_mpi_type), "MPI_Type_commit");
		send_mpi_types_.push_back(new_mpi_type);

		int type_size = 0;
		check_error(MPI_Type_size(new_mpi_type, &type_size), "MPI_Type_size");
		send_mpi_type_sizes_.push_back(static_cast<size_t>(type_size));

		DEBUG_ONLY(
		    std::cout << "communicator::communicator(..): created MPI_Type_indexed to send from rank " << rank_ << " to " << destination
		              << " of "
		              << boost::format("%.2f MiB") % (type_size / 1024.0 / 1024.0)
		              << " in " << block_lengths.size() << " blocks" << std::endl;
		)
		DEBUG_ONLY(
			for (size_t i = 0; i < displacements.size(); ++i) {
				std::cout << "communicator::communicator(..): rank " << rank_ << " to " << destination << ": block "
				          << i << ": l = " << block_lengths[i] << ", d = " << displacements[i] << std::endl;
				assert(((displacements[i] + block_lengths[i])) <= (partition.compute_nodes()));
			}
		)
		ASSERT_ONLY( last_destination = destination; )
	}
	assert(destinations_.size() == send_counts_.size() && send_counts_.size() == send_displacements_.size() && send_displacements_.size() == send_mpi_types_.size());


	// setup recv/in/sources stuff (simple, contiguous sets of nodes needed by this rank sent from other ranks, matched by a non-contiguous send data type on another rank)

	ASSERT_ONLY( rank_t last_source = -1; )
	int node_offset = 0;
	for (const auto& source : sources_)
	{
		assert(source > last_source); // make sure the expected rank order is valid, i.e. numerical order

		recv_counts_.push_back(1); // NOTE: currently exactly one instance per type is required
		const int node_count = partition.rank_to_in_neighbours().at(source).size();
		const int displacement_byte = node_offset * inst.size_hierarchy_top_byte();
		recv_displacements_.push_back(displacement_byte);
		node_offset += node_count;

		// data is already aligned contiguously but might vary in size for each neighbour partition
		// contiguous buffer containing the amount of hierarchy nodes received by a certain neighbour rank
		MPI_Datatype new_mpi_type;
		check_error(MPI_Type_contiguous(node_count, hierarchy_node_mpi_type_, &new_mpi_type), "MPI_Type_contiguous");
		check_error(MPI_Type_commit(&new_mpi_type), "MPI_Type_commit");
		recv_mpi_types_.push_back(new_mpi_type);

		int type_size = 0;
		check_error(MPI_Type_size(new_mpi_type, &type_size), "MPI_Type_size");
		recv_mpi_type_sizes_.push_back(static_cast<size_t>(type_size));

		DEBUG_ONLY(
			std::cout << "communicator::communicator(..): created MPI_Type_contiguous to recv on rank " << rank_ << " from " << source
			          << " of "
			          << boost::format("%.2f MiB") % (type_size / 1024.0 / 1024.0)
			          << " with " << node_count << " hierarchy nodes, displaced by " << displacement_byte << " byte" << std::endl;
		)

		ASSERT_ONLY( last_source = source; )
	}
	assert(sources_.size() == recv_counts_.size() && recv_counts_.size() == recv_displacements_.size() && recv_displacements_.size() == recv_mpi_types_.size());

	// exchange data type sizes

	// prepare vector with this rank's values
	std::vector<size_t> send_row(ranks_, 0); // row of send sizes with index = rank for this rank
	for (size_t i = 0; i < destinations_.size(); ++i)
		send_row[destinations_[i]] = send_mpi_type_sizes_[i];

	// exchange rows between all nodes
	const size_t data_size = ranks_ * sizeof(size_t); // one row of send_size_matrix_
	check_error(MPI_Allgather(send_row.data(), //const void *sendbuf
	                          data_size, //int sendcount
	                          MPI_BYTE, //MPI_Datatype sendtype
	                          send_size_matrix_.data(), //void *recvbuf
	                          data_size, //int recvcount
		                      MPI_BYTE,
	                          MPI_COMM_WORLD), // MPI_Comm comm
	            "MPI_Allgather");
}

communicator::~communicator()
{
	check_error(MPI_Type_free(&hierarchy_node_mpi_type_), "MPI_Type_free");

	for (auto& type : send_mpi_types_)
		check_error(MPI_Type_free(&type), "MPI_Type_free");

	for (auto& type : recv_mpi_types_)
		check_error(MPI_Type_free(&type), "MPI_Type_free");
}

void communicator::barrier()
{
	check_error(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier");
}

void communicator::neighbour_exchange(void* send_buffer, void* recv_buffer)
{
	// actual data transfer
	DEBUG_ONLY( std::cout << "communicator::neighbour_exchange(): calling MPI_Neighbor_alltoallw(..) on rank " << rank_ << std::endl; )
	check_error(MPI_Neighbor_alltoallw(send_buffer, // const void *sendbuf,
	                                   send_counts_.data(), // const int sendcounts[],
	                                   send_displacements_.data(), // send_displacements_.data(), // const MPI_Aint sdispls[],
	                                   send_mpi_types_.data(), // const MPI_Datatype sendtypes[],
	                                   recv_buffer, // void *recvbuf,
	                                   recv_counts_.data(), // const int recvcounts[],
	                                   recv_displacements_.data(), // const MPI_Aint rdispls[],
	                                   recv_mpi_types_.data(), // const MPI_Datatype recvtypes[],
	                                   neighbour_comm_), // MPI_Comm comm)
	            "MPI_Neighbor_alltoallw");
}

void communicator::write_neighbour_exchange_table(std::ostream& os)
{
	// NOTE: we produce human readable output here, use the _raw version for further processing
	//const char delimiter = default_delimiter;
	const char delimiter = ' ';

	// formatting
	std::stringstream entry_format;
	entry_format << "%9.2f MiB";  // 9-3 digits before the dot
	std::stringstream entry_sample;
	entry_sample << boost::format(entry_format.str()) % 0.0; // to get the length
	std::stringstream first_row_format; // format string for rank header
	first_row_format << "%" << entry_sample.str().length() << "d"; // header spacing must match entry_format
	std::string first_col_head = "send_rank";
	std::string last_row_begin = " recv_sum"; // same length as above
	std::string last_col_head = "send_sum";
	last_col_head = std::string(entry_sample.str().length() - last_col_head.length(), ' ') + last_col_head; // prepend spaces to match entry
	std::stringstream first_col_format; // format string
	first_col_format << "%" << first_col_head.length() << "d";
	std::stringstream last_col_format; // format string
	last_col_format << "%" << last_col_head.length() << "d";

	auto size_to_mib = [&](size_t size) {
		std::stringstream ss;
		ss << (boost::format(entry_format.str()) % (size / 1024.0 / 1024.0));
		return ss.str();
	};

	// header with ranks
	os << first_col_head;
	for (rank_t r = 0; r < ranks_; ++r)
		os << delimiter << boost::format(first_row_format.str()) % r;
	os << delimiter << last_col_head << "\n";

	// rows
	std::vector<size_t> col_sums(ranks_, 0);
	for (rank_t r = 0; r < ranks_; ++r) { // rows
		os << boost::format(first_col_format.str()) % r; // current rank in first column
		size_t row_sum = 0;
		for (rank_t c = 0; c < ranks_; ++c) { // columns in each row
			size_t& current_size = send_size_matrix_.at(r, c);
			os << delimiter << size_to_mib(current_size); // write value
			row_sum += current_size; // count row sum of current row
			col_sums[c] += current_size; // count column sums
		}
		os << delimiter << size_to_mib(row_sum) << "\n"; // last column contains row sum
	}

	// last row contains column sums, i.e. recv counts, for rank in header of that column
	os << last_row_begin;
	size_t data_volume = 0;
	for (rank_t r = 0; r < ranks_; ++r) { // columns of last row
		size_t current_col_sum = col_sums[r];
		os << delimiter << size_to_mib(current_col_sum);
		data_volume += current_col_sum;
	}
	os << delimiter << size_to_mib(data_volume) << std::endl;
}

void communicator::write_neighbour_exchange_table_raw(std::ostream& os)
{
	send_size_matrix_.print_raw(os, default_delimiter);
}


} // namespace heom
