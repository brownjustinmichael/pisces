/*!**********************************************************************
 * \file messenger.cpp
 * /Users/justinbrown/Dropbox/spectral_element
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
	/*
		TODO This module could use a lot of cleaning...
	*/
	
	template <class datatype>
	messenger <datatype>::messenger (int* argc, char*** argv, int n_boundaries) {
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
	
	template <class datatype>
	messenger <datatype>::~messenger () {
#ifdef _MPI
		MPI::Finalize ();
#endif // _MPI
	}
	
	template <class datatype>
	void messenger <datatype>::add_boundary (int edge, int process) {
		process_queue [edge_to_index (send_mode, edge)] = process;
		process_queue [edge_to_index (recv_mode, edge)] = process;
	}
	
	template <class datatype>
	void messenger <datatype>::send (int n, datatype* data, int edge) {
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
		data_check ();
	}
	
	template <class datatype>
	void messenger <datatype>::recv (int n, datatype* data, int edge) {
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
		data_check ();
	}

	template <>
	void messenger <double>::data_check () {
		TRACE ("Running messenger check.");
		while (data_queue [data_iter]) {
			if (process_queue [data_iter] == -1 || n_queue [data_iter] == 0) {
				if (index_to_mode (data_iter) == recv_mode) {
					utils::scale (n_queue [data_iter], 0.0, data_queue [data_iter]);
				}
			} else {
				if (index_to_mode (data_iter) == send_mode) {
#ifdef _MPI
					TRACE ("Performing send.");
					MPI::COMM_WORLD.Send (data_queue [data_iter], n_queue [data_iter], MPI::DOUBLE, process_queue [data_iter], data_iter);
#else // _MPI
					FATAL ("Send used without MPI environment. Exiting.");
					throw 0;
#endif // _MPI
				} else {
#ifdef _MPI
					TRACE ("Performing recv.");
					MPI::COMM_WORLD.Recv (data_queue [data_iter], n_queue [data_iter], MPI::DOUBLE, process_queue [data_iter], data_iter);
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
	
	template <>
	void messenger <float>::data_check () {
		TRACE ("Running messenger check.");
		while (data_queue [data_iter]) {
			if (process_queue [data_iter] == -1 || n_queue [data_iter] == 0) {
				if (index_to_mode (data_iter) == recv_mode) {
					utils::scale (n_queue [data_iter], 0.0, data_queue [data_iter]);
				}
			} else {
				if (index_to_mode (data_iter) == send_mode) {
#ifdef _MPI
					TRACE ("Performing send.");
					MPI::COMM_WORLD.Send (data_queue [data_iter], n_queue [data_iter], MPI::FLOAT, process_queue [data_iter], data_iter);
#else // _MPI
					FATAL ("Send used without MPI environment. Exiting.");
					throw 0;
#endif // _MPI
				} else {
#ifdef _MPI
					TRACE ("Performing recv.");
					MPI::COMM_WORLD.Recv (data_queue [data_iter], n_queue [data_iter], MPI::FLOAT, process_queue [data_iter], data_iter);
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
	void messenger <datatype>::send (int n, int* data, int edge) {
		TRACE ("Adding message to queue.");
		if (int_data_queue [edge_to_index (send_mode, edge)]) {
			FATAL ("Message already in queue at this edge.");
			throw 0;
		}
		if (!data) {
			if (n == 0) {
				data = (int*) 0x01;
			} else {
				WARN ("Nonzero null pointer specified for messenger operation");
			}
		}
		int_data_queue [edge_to_index (send_mode, edge)] = data;
		n_queue [edge_to_index (send_mode, edge)] = n;
		int_check ();
	}
	
	template <class datatype>
	void messenger <datatype>::recv (int n, int* data, int edge) {
		TRACE ("Adding message to queue.");
		if (int_data_queue [edge_to_index (recv_mode, edge)]) {
			FATAL ("Message already in queue at this edge.");
			throw 0;
		}
		if (!data) {
			if (n == 0) {
				data = (int*) 0x01;
			} else {
				WARN ("Nonzero null pointer specified for messenger operation");
			}
		}
		int_data_queue [edge_to_index (recv_mode, edge)] = data;
		n_queue [edge_to_index (recv_mode, edge)] = n;
		int_check ();
	}
	
	template <class datatype>
	void messenger <datatype>::int_check () {
		TRACE ("Running messenger check.");
		while (int_data_queue [int_iter]) {
			if (process_queue [int_iter] == -1 || n_queue [int_iter] == 0) {
				if (index_to_mode (int_iter) == recv_mode) {
					TRACE ("Zeroing recv array.");
					for (int i = 0; i < n_queue [int_iter]; ++i) {
						int_data_queue [int_iter] [i] = 0;
					}
				}
			} else {
				if (index_to_mode (int_iter) == send_mode) {
#ifdef _MPI
					TRACE ("Performing send.");
					MPI::COMM_WORLD.Send (int_data_queue [int_iter], n_queue [int_iter], MPI::INT, process_queue [int_iter], int_iter);
#else // _MPI
					FATAL ("Send used without MPI environment. Exiting.");
					throw 0;
#endif // _MPI
				} else {
#ifdef _MPI
					TRACE ("Performing recv.");
					MPI::COMM_WORLD.Recv (int_data_queue [int_iter], n_queue [int_iter], MPI::INT, process_queue [int_iter], int_iter);
#else // _MPI
					FATAL ("Recv used without MPI environment. Exiting.");
					throw 0;
#endif // _MPI
				}
			}
			int_data_queue [int_iter] = NULL;
			if (int_iter + 1 >= (int) int_data_queue.size ()) {
				int_iter = 0;
			} else {
				++int_iter;
			}
		}
		TRACE ("Check complete.");
	}
	
	template <>
	void messenger <float>::min (float* data) {
		if (np != 1) {
			if (id == 0 && np > (int) buffer.size ()) {
				buffer.resize (np);
			}
#ifdef _MPI
			MPI::COMM_WORLD.Gather (data, 1, MPI::FLOAT, &buffer [0], 1, MPI::DOUBLE, 0);
#else // _MPI
			FATAL ("Gather used without MPI environment. Exiting.");
			throw 0;
#endif // _MPI
			if (id == 0) {
				*data = *std::min_element (buffer.begin (), buffer.end ());
			}
#ifdef _MPI
			MPI::COMM_WORLD.Bcast (data, 1, MPI::FLOAT, 0);
#else // _MPI
			FATAL ("Comm_world used without MPI environment. Exiting.");
			throw 0;
#endif // _MPI
		}
	}
	
	template <>
	void messenger <double>::min (double* data) {
		if (np != 1) {
			if (id == 0 && np > (int) buffer.size ()) {
				buffer.resize (np);
			}
#ifdef _MPI
			MPI::COMM_WORLD.Gather (data, 1, MPI::DOUBLE, &buffer [0], 1, MPI::DOUBLE, 0);
#else // _MPI
			FATAL ("Gather used without MPI environment. Exiting.");
			throw 0;
#endif // _MPI
			if (id == 0) {
				*data = *std::min_element (buffer.begin (), buffer.end ());
			}
#ifdef _MPI
			MPI::COMM_WORLD.Bcast (data, 1, MPI::DOUBLE, 0);
#else // _MPI
			FATAL ("Comm_world used without MPI environment. Exiting.");
			throw 0;
#endif // _MPI
		}
	}
	
	template <class datatype>
	bool messenger <datatype>::bool_and (bool boolean) {
		int temp;
		if (boolean) {
			temp = 1;
		} else {
			temp = 0;
		}
		if (np != 1) {
			if (id == 0 && np > (int) int_buffer.size ()) {
				int_buffer.resize (np);
			}
#ifdef _MPI
			MPI::COMM_WORLD.Gather (&temp, 1, MPI::INT, &int_buffer [0], 1, MPI::INT, 0);
#else // _MPI
			FATAL ("Gather used without MPI environment. Exiting.");
			throw 0;
#endif // _MPI
			if (id == 0) {
				for (int i = 0; i < np; i++) {
					if (!int_buffer [i]) {
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
	
	template class messenger <double>;
	template class messenger <float>;	
} /* bases */
