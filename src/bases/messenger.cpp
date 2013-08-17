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
#include "mpi.h"
#include "../config.hpp"
#include "../utils/utils.hpp"

namespace bases
{	
	messenger::messenger (int* argc, char*** argv, int n_boundaries) {
		MPI::Init (*argc, *argv);
		np = MPI::COMM_WORLD.Get_size ();
		id = MPI::COMM_WORLD.Get_rank ();
		assert (n_boundaries % 2 == 0);
		double_data_queue.resize (2 * n_boundaries);
		int_data_queue.resize (2 * n_boundaries);
		n_queue.resize (2 * n_boundaries);
		process_queue.resize (2 * n_boundaries, -1);
	}
	
	messenger::~messenger () {
		MPI::Finalize ();
	}
	
	void messenger::add_boundary (int edge, int process) {
		process_queue [edge_to_index (send_mode, edge)] = process;
		process_queue [edge_to_index (recv_mode, edge)] = process;
	}
	
	void messenger::send (int n, double* data, int edge) {
		if (double_data_queue [edge_to_index (send_mode, edge)]) {
			FATAL ("Message already in queue at this edge.");
			throw 0;
		}
		if (!data) {
			if (n == 0) {
				data = (double*) 0x01;
			} else {
				WARN ("Nonzero null pointer specified for messenger operation");
			}
		}
		double_data_queue [edge_to_index (send_mode, edge)] = data;
		n_queue [edge_to_index (send_mode, edge)] = n;
		double_check ();
	}
	
	void messenger::recv (int n, double* data, int edge) {
		if (double_data_queue [edge_to_index (recv_mode, edge)]) {
			FATAL ("Message already in queue at this edge.");
			throw 0;
		}
		if (!data) {
			if (n == 0) {
				data = (double*) 0x01;
			} else {
				WARN ("Nonzero null pointer specified for messenger operation");
			}
		}
		double_data_queue [edge_to_index (recv_mode, edge)] = data;
		n_queue [edge_to_index (recv_mode, edge)] = n;
		double_check ();
	}

	void messenger::double_check () {
		while (double_data_queue [double_iter]) {
			if (process_queue [double_iter] == -1 || n_queue [double_iter] == 0) {
				if (index_to_mode (double_iter) == recv_mode) {
					utils::scale (n_queue [double_iter], 0.0, double_data_queue [double_iter]);
				}
			} else {
				if (index_to_mode (double_iter) == send_mode) {
					MPI::COMM_WORLD.Send (double_data_queue [double_iter], n_queue [double_iter], MPI::DOUBLE, process_queue [double_iter], double_iter);
				} else {
					MPI::COMM_WORLD.Recv (double_data_queue [double_iter], n_queue [double_iter], MPI::DOUBLE, process_queue [double_iter], double_iter);
				}
			}
			double_data_queue [double_iter] = NULL;
			if (double_iter + 1 >= (int) double_data_queue.size ()) {
				double_iter = 0;
			} else {
				++double_iter;
			}
		}
	}

	void messenger::send (int n, int* data, int edge) {
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
	
	void messenger::recv (int n, int* data, int edge) {
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
	
	void messenger::int_check () {
		while (int_data_queue [int_iter]) {
			if (process_queue [int_iter] == -1 || n_queue [int_iter] == 0) {
				if (index_to_mode (int_iter) == recv_mode) {
					for (int i = 0; i < n_queue [int_iter]; ++i) {
						int_data_queue [int_iter] [i] = 0;
					}
				}
			} else {
				if (index_to_mode (int_iter) == send_mode) {
					MPI::COMM_WORLD.Send (int_data_queue [int_iter], n_queue [int_iter], MPI::INT, process_queue [int_iter], int_iter);
				} else {
					MPI::COMM_WORLD.Recv (int_data_queue [int_iter], n_queue [int_iter], MPI::INT, process_queue [int_iter], int_iter);
				}
			}
			int_data_queue [int_iter] = NULL;
			if (int_iter + 1 >= (int) int_data_queue.size ()) {
				int_iter = 0;
			} else {
				++int_iter;
			}
		}
	}
	
	void messenger::min (double* data) {
		if (np != 1) {
			if (id == 0 && np > (int) buffer.size ()) {
				buffer.resize (np);
			}
			MPI::COMM_WORLD.Gather (data, 1, MPI::DOUBLE, &buffer [0], 1, MPI::DOUBLE, 0);
			if (id == 0) {
				*data = *std::min_element (buffer.begin (), buffer.end ());
			}
			MPI::COMM_WORLD.Bcast (data, 1, MPI::DOUBLE, 0);
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
			if (id == 0 && np > (int) int_buffer.size ()) {
				int_buffer.resize (np);
			}
			MPI::COMM_WORLD.Gather (&temp, 1, MPI::INT, &int_buffer [0], 1, MPI::INT, 0);
			if (id == 0) {
				for (int i = 0; i < np; i++) {
					if (!int_buffer [i]) {
						temp = 0;
						break;
					}
				}
				MPI::COMM_WORLD.Bcast (&temp, 1, MPI::INT, 0);
			}
		}
		if (temp) {
			return true;
		} else {
			return false;
		}
	}
		
} /* bases */
