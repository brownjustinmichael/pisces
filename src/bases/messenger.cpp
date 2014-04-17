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
		MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
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
		stati.resize (np);
	}
	
	messenger::~messenger () {
		// printf ("Destroying bases messenger\n");
		kill_all ();
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
	void messenger::send (int n, datatype *data, int process, int tag) {
#ifdef _MPI
		MPI::COMM_WORLD.Send (data, n, mpi_type <datatype> (), process, tag);
#else // _MPI
		FATAL ("Send used without MPI environment. Exiting.");
		// throw 0;
#endif // _MPI
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
	void messenger::recv (int n, datatype *data, int process, int tag) {
		/*
			TODO Make check_two
		*/
#ifdef _MPI
		MPI::COMM_WORLD.Recv (data, n, mpi_type <datatype> (), process, tag);
#else // _MPI
		FATAL ("Send used without MPI environment. Exiting.");
		// throw 0;
#endif // _MPI
	}

	template <class datatype>
	void messenger::data_check () {
		int flags = mpi_all_clear;
		check_all (&flags);
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
					// throw 0;
#endif // _MPI
				} else {
#ifdef _MPI
					TRACE ("Performing recv.");
					MPI::COMM_WORLD.Recv (data_queue [data_iter], n_queue [data_iter], mpi_type <datatype> (), process_queue [data_iter], data_iter);
#else // _MPI
					FATAL ("Send used without MPI environment. Exiting.");
					// throw 0;
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
		int flags = mpi_all_clear;
		check_all (&flags);
		if (np != 1) {
			std::vector <datatype> buffer (np);
#ifdef _MPI
			MPI::COMM_WORLD.Gather (data, 1, mpi_type <datatype> (), &buffer [0], 1, mpi_type <datatype> (), 0);
#else // _MPI
			FATAL ("Gather used without MPI environment. Exiting.");
			// throw 0;
#endif // _MPI
			if (id == 0) {
				*data = *std::min_element (buffer.begin (), buffer.end ());
			}
#ifdef _MPI
			MPI::COMM_WORLD.Bcast (data, 1, mpi_type <datatype> (), 0);
#else // _MPI
			FATAL ("Comm_world used without MPI environment. Exiting.");
			// throw 0;
#endif // _MPI
		}
	}
	
	void messenger::check_all (int *flags) {
		TRACE ("Checking...")
#ifdef _MPI
		MPI::COMM_WORLD.Gather (flags, 1, mpi_type <int> (), &stati [0], 1, mpi_type <int> (), 0);
#else // _MPI
		FATAL ("Gather used without MPI environment. Exiting.");
		// throw 0;
#endif // _MPI
		if (id == 0) {
			for (int i = 0; i < np; ++i) {
				*flags |= stati [i];
			}
		}
#ifdef _MPI
		MPI::COMM_WORLD.Bcast (flags, 1, mpi_type <int> (), 0);
#else // _MPI
		FATAL ("Comm_world used without MPI environment. Exiting.");
		// throw 0;
#endif // _MPI
		if (*flags & mpi_fatal) {
			throw 0;
		}
		TRACE ("Check complete.");
	}
	
	void messenger::kill_all () {
		int flags = mpi_fatal;
#ifdef _MPI
		MPI::COMM_WORLD.Bcast (&flags, 1, mpi_type <int> (), 0);
#else // _MPI
		FATAL ("Comm_world used without MPI environment. Exiting.");
		// throw 0;
#endif // _MPI
	}
	
	template <class datatype>
	void messenger::gather (int n, datatype* data_in, datatype* data_out) {
		int flags = mpi_all_clear;
		check_all (&flags);
		if (!data_out) {
			data_out = data_in;
		}
#ifdef _MPI
		if (id == 0 && data_out == data_in) {
			MPI::COMM_WORLD.Gather (MPI_IN_PLACE, n, mpi_type <datatype> (), data_in, n, mpi_type <datatype> (), 0);
		} else {
			MPI::COMM_WORLD.Gather (data_in, n, mpi_type <datatype> (), data_out, n, mpi_type <datatype> (), 0);
		}
#else // _MPI
		FATAL ("Gather used without MPI environment. Exiting.");
		// throw 0;
#endif // _MPI
	}
	
	template <class datatype>
	void messenger::bcast (int n, datatype* data_in) {
		int flags = mpi_all_clear;
		check_all (&flags);
#ifdef _MPI
		if (flags == mpi_all_clear) {
			MPI::COMM_WORLD.Bcast (data_in, n, mpi_type <datatype> (), 0);
		}
		/*
			TODO Use something like this instead of assuming check_all will always succeed or throw an exception?
		*/
#else // _MPI
		FATAL ("Gather used without MPI environment. Exiting.");
		// throw 0;
#endif // _MPI
	}
	
	template <class datatype>
	void messenger::allgather (int n, datatype* data_in, datatype* data_out) {
		int flags = mpi_all_clear;
		check_all (&flags);
		if (!data_out) {
			data_out = data_in;
		}
#ifdef _MPI
		if (data_out == data_in) {
			MPI::COMM_WORLD.Allgather (MPI_IN_PLACE, n, mpi_type <datatype> (), data_in, n, mpi_type <datatype> ());
		} else {
			MPI::COMM_WORLD.Allgather (data_in, n, mpi_type <datatype> (), data_out, n, mpi_type <datatype> ());
		}
#else // _MPI
		FATAL ("Gather used without MPI environment. Exiting.");
		// throw 0;
#endif // _MPI
	}
	
	template <class datatype>
	void messenger::gatherv (int n, datatype* data_in, int *ns, datatype* data_out) {
		int flags = mpi_all_clear;
		check_all (&flags);
		if (!data_out) {
			data_out = data_in;
		}
#ifdef _MPI
		std::vector <int> displs;
		if (id == 0) {
			displs.resize (np);
			for (int i = 1; i < np; ++i) {
				displs [i] = displs [i - 1] + ns [i - 1];
			}
		}
		if (id == 0 && data_out == data_in) {
			MPI::COMM_WORLD.Gatherv (MPI_IN_PLACE, n, mpi_type <datatype> (), data_in, ns, &displs [0], mpi_type <datatype> (), 0);
		} else {
			MPI::COMM_WORLD.Gatherv (data_in, n, mpi_type <datatype> (), data_out, ns, &displs [0], mpi_type <datatype> (), 0);
		}
#else // _MPI
		FATAL ("Gather used without MPI environment. Exiting.");
		// throw 0;
#endif // _MPI
	}
	
	template <class datatype>
	void messenger::allgatherv (int n, datatype* data_in, int *ns, datatype* data_out) {
		int flags = mpi_all_clear;
		check_all (&flags);
		if (!data_out) {
			data_out = data_in;
		}
#ifdef _MPI
		std::vector <int> displs;
		displs.resize (np);
		for (int i = 1; i < np; ++i) {
			displs [i] = displs [i - 1] + ns [i - 1];
		}
		if (data_out == data_in) {
			MPI::COMM_WORLD.Allgatherv (MPI_IN_PLACE, n, mpi_type <datatype> (), data_in, ns, &displs [0], mpi_type <datatype> ());
		} else {
			MPI::COMM_WORLD.Allgatherv ((void *) data_in, n, mpi_type <datatype> (), (void *) data_out, ns, &displs [0], mpi_type <datatype> ());
		}
#else // _MPI
		FATAL ("Gather used without MPI environment. Exiting.");
// 		// throw 0;
#endif // _MPI
	}
	
	template <class datatype>
	void messenger::scatter (int n, datatype* data_in, datatype* data_out) {
		int flags = mpi_all_clear;
		check_all (&flags);
		if (!data_out) {
			data_out = data_in;
		}
#ifdef _MPI
		if (id == 0 && data_out == data_in) {
			MPI::COMM_WORLD.Scatter (data_in, n, mpi_type <datatype> (), MPI_IN_PLACE, n, mpi_type <datatype> (), 0);
		} else {
			MPI::COMM_WORLD.Scatter (data_in, n, mpi_type <datatype> (), data_out, n, mpi_type <datatype> (), 0);
		}
#else // _MPI
		FATAL ("Scatter used without MPI environment. Exiting.");
		// throw 0;
#endif // _MPI
	}
	
	template <class datatype>
	void messenger::scatterv (int n, datatype* data_in, int* ns, datatype* data_out) {
		int flags = mpi_all_clear;
		check_all (&flags);
		if (!data_out) {
			data_out = data_in;
		}
#ifdef _MPI
		std::vector <int> displs;
		if (id == 0) {
			displs.resize (np);
			for (int i = 1; i < np; ++i) {
				displs [i] = displs [i - 1] + ns [i - 1];
			}
		}
		if (id == 0 && data_out == data_in) {
			MPI::COMM_WORLD.Scatterv (data_in, ns, &displs [0], mpi_type <datatype> (), MPI_IN_PLACE, n, mpi_type <datatype> (), 0);
		} else {
			MPI::COMM_WORLD.Scatterv (data_in, ns, &displs [0], mpi_type <datatype> (), data_out, n, mpi_type <datatype> (), 0);
		}
#else // _MPI
		FATAL ("Scatter used without MPI environment. Exiting.");
		// throw 0;
#endif // _MPI
	}
	
	bool messenger::bool_and (bool boolean) {
		int flags = mpi_all_clear;
		check_all (&flags);
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
			// throw 0;
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
				// throw 0;
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
	
	template void messenger::send <double> (int n, double* data, int process, int tag);
	template void messenger::send <float> (int n, float* data, int process, int tag);
	template void messenger::send <int> (int n, int* data, int process, int tag);
	
	template void messenger::recv <double> (int n, double* data, int process, int tag);
	template void messenger::recv <float> (int n, float* data, int process, int tag);
	template void messenger::recv <int> (int n, int* data, int process, int tag);
	
	template void messenger::recv <double> (int n, double* data, int edge);
	template void messenger::recv <float> (int n, float* data, int edge);
	template void messenger::recv <int> (int n, int* data, int edge);
	
	template void messenger::data_check <double> ();
	template void messenger::data_check <float> ();
	template void messenger::data_check <int> ();
	
	template void messenger::min <double> (double* data);
	template void messenger::min <float> (float* data);
	template void messenger::min <int> (int* data);

	template void messenger::gather <double> (int n, double* data_in, double* data_out);
	template void messenger::gather <float> (int n, float* data_in, float* data_out);
	template void messenger::gather <int> (int n, int* data_in, int* data_out);
	
	template void messenger::bcast <double> (int n, double* data_in);
	template void messenger::bcast <float> (int n, float* data_in);
	template void messenger::bcast <int> (int n, int* data_in);
	
	template void messenger::allgather <double> (int n, double* data_in, double* data_out);
	template void messenger::allgather <float> (int n, float* data_in, float* data_out);
	template void messenger::allgather <int> (int n, int* data_in, int* data_out);
	
	template void messenger::gatherv <double> (int n, double* data_in, int *ns, double* data_out);
	template void messenger::gatherv <float> (int n, float* data_in, int *ns, float* data_out);
	template void messenger::gatherv <int> (int n, int* data_in, int *ns, int* data_out);
	
	template void messenger::allgatherv <double> (int n, double* data_in, int *ns, double* data_out);
	template void messenger::allgatherv <float> (int n, float* data_in, int *ns, float* data_out);
	template void messenger::allgatherv <int> (int n, int* data_in, int *ns, int* data_out);
	
	template void messenger::scatter <double> (int n, double* data_in, double* data_out);
	template void messenger::scatter <float> (int n, float* data_in, float* data_out);
	template void messenger::scatter <int> (int n, int* data_in, int* data_out);
	
	template void messenger::scatterv <double> (int n, double* data_in, int *ns, double* data_out);
	template void messenger::scatterv <float> (int n, float* data_in, int *ns, float* data_out);
	template void messenger::scatterv <int> (int n, int* data_in, int *ns, int* data_out);
} /* bases */
