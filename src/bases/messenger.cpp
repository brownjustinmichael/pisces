/*!**********************************************************************
 * \file messenger.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-08-17.
* Copyright 2013 Justin Brown. All rights reserved.
************************************************************************/

#include "messenger.hpp"
#include <cassert>
#include <algorithm>
#include <vector>
#ifdef _MPI
#include "mpi.h"
#endif // _MPI
#include "../config.hpp"
#include "../utils/utils.hpp"

namespace bases
{
#ifdef _MPI
	template <class datatype>
	MPI_Datatype mpi_type ();
	
	template <>
	MPI_Datatype mpi_type <double> () {
		return MPI::DOUBLE;
	}
	
	template <>
	MPI_Datatype mpi_type <float> () {
		return MPI::REAL;
	}
	
	template <>
	MPI_Datatype mpi_type <int> () {
		return MPI::INT;
	}
#endif // _MPI
	
	messenger::messenger (int* argc, char*** argv, int n_boundaries) {
#ifdef _MPI
		MPI::Init (*argc, *argv);
		np = MPI::COMM_WORLD.Get_size ();
		id = MPI::COMM_WORLD.Get_rank ();
#else
		np = 1;
		id = 0;
#endif // _MPI
		assert (n_boundaries % 2 == 0);
		data_iter = 0;
		int_iter = 0;
		data_queue.resize (2 * n_boundaries);
		int_data_queue.resize (2 * n_boundaries);
		n_queue.resize (2 * n_boundaries);
		process_queue.resize (2 * n_boundaries, -1);
	}
	
	messenger::~messenger () {
#ifdef _MPI
		MPI::Finalize ();
#endif // _MPI
	}
	
	void messenger::add_boundary (int edge, int process) {
		process_queue [edge_to_index (send_mode, edge)] = process;
		process_queue [edge_to_index (recv_mode, edge)] = process;
	}
	
	template <class datatype>
	void messenger::send (int n, datatype* data, int edge) {
		TRACE ("Adding message to queue.");
		if (data_queue [edge_to_index (send_mode, edge)]) {
			FATAL ("Message already in queue at this edge.");
			throw 0;
		}
		if (!data) {
			if (n == 0) {
				data = (datatype*) 0x01;
			} else {
				WARN ("Nonzero null pointer specified for messenger operation");
			}
		}
		data_queue [edge_to_index (send_mode, edge)] = data;
		n_queue [edge_to_index (send_mode, edge)] = n;
		data_check <datatype> ();
	}
	
	template <class datatype>
	void messenger::recv (int n, datatype* data, int edge) {
		TRACE ("Adding message to queue.");
		if (data_queue [edge_to_index (recv_mode, edge)]) {
			FATAL ("Message already in queue at this edge.");
			throw 0;
		}
		if (!data) {
			if (n == 0) {
				data = (datatype*) 0x01;
			} else {
				WARN ("Nonzero null pointer specified for messenger operation");
			}
		}
		data_queue [edge_to_index (recv_mode, edge)] = data;
		n_queue [edge_to_index (recv_mode, edge)] = n;
		data_check <datatype> ();
	}

	template <class datatype>
	void messenger::data_check () {
		TRACE ("Running messenger check.");
		while (data_queue [data_iter]) {
			if (process_queue [data_iter] == -1 || n_queue [data_iter] == 0) {
				if (index_to_mode (data_iter) == recv_mode) {
					utils::scale (n_queue [data_iter], (datatype) 0, (datatype*) data_queue [data_iter]);
				}
			} else {
				if (index_to_mode (data_iter) == send_mode) {
#ifdef _MPI
					TRACE ("Performing send.");
					MPI::COMM_WORLD.Send (data_queue [data_iter], n_queue [data_iter], mpi_type <datatype> (), process_queue [data_iter], data_iter);
#else // _MPI
					FATAL ("Send used without MPI environment. Exiting.");
					throw 0;
#endif // _MPI
				} else {
#ifdef _MPI
					TRACE ("Performing recv.");
					MPI::COMM_WORLD.Recv (data_queue [data_iter], n_queue [data_iter], mpi_type <datatype> (), process_queue [data_iter], data_iter);
#else // _MPI
					FATAL ("Send used without MPI environment. Exiting.");
					throw 0;
#endif // _MPI
				}
			}
			data_queue [data_iter] = NULL;
			if (data_iter + 1 >= (int) data_queue.size ()) {
				data_iter = 0;
			} else {
				++data_iter;
			}
		}
		TRACE ("Check complete.");
	}
	
	template <class datatype>
	void messenger::min (datatype* data) {
		if (np != 1) {
			std::vector <datatype> buffer (np);
#ifdef _MPI
			MPI::COMM_WORLD.Gather (data, 1, mpi_type <datatype> (), &buffer [0], 1, mpi_type <datatype> (), 0);
#else // _MPI
			FATAL ("Gather used without MPI environment. Exiting.");
			throw 0;
#endif // _MPI
			if (id == 0) {
				*data = *std::min_element (buffer.begin (), buffer.end ());
			}
#ifdef _MPI
			MPI::COMM_WORLD.Bcast (data, 1, mpi_type <datatype> (), 0);
#else // _MPI
			FATAL ("Comm_world used without MPI environment. Exiting.");
			throw 0;
#endif // _MPI
		}
	}
	
	bool messenger::bool_and (bool boolean) {
		int temp;
		if (boolean) {
			temp = 1;
		} else {
			temp = 0;
		}
		if (np != 1) {
			std::vector <int> buffer;
			if (id == 0 && np > (int) buffer.size ()) {
				buffer.resize (np);
			}
#ifdef _MPI
			MPI::COMM_WORLD.Gather (&temp, 1, MPI::INT, &buffer [0], 1, MPI::INT, 0);
#else // _MPI
			FATAL ("Gather used without MPI environment. Exiting.");
			throw 0;
#endif // _MPI
			if (id == 0) {
				for (int i = 0; i < np; i++) {
					if (!buffer [i]) {
						temp = 0;
						break;
					}
				}
#ifdef _MPI
				MPI::COMM_WORLD.Bcast (&temp, 1, MPI::INT, 0);
#else // _MPI
				FATAL ("Bcast used without MPI environment. Exiting.");
				throw 0;
#endif // _MPI
			}
		}
		if (temp) {
			return true;
		} else {
			return false;
		}
	}
	
	template void messenger::send <double> (int n, double* data, int edge);
	template void messenger::send <float> (int n, float* data, int edge);
	template void messenger::send <int> (int n, int* data, int edge);
	
	template void messenger::recv <double> (int n, double* data, int edge);
	template void messenger::recv <float> (int n, float* data, int edge);
	template void messenger::recv <int> (int n, int* data, int edge);
	
	template void messenger::data_check <double> ();
	template void messenger::data_check <float> ();
	template void messenger::data_check <int> ();
	
	template void messenger::min <double> (double* data);
	template void messenger::min <float> (float* data);
	template void messenger::min <int> (int* data);
} /* bases */
