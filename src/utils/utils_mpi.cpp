/*!***********************************************************************
 * \file bases/message.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "messenger.hpp"
#include "../config.hpp"
#include <algorithm>
#include "mpi.h"
#include "utils.hpp"

#define OPENMPI

namespace utils
{
	messenger::messenger (int* argc, char*** argv) {
		MPI::Init (*argc, *argv);
		np = MPI::COMM_WORLD.Get_size ();
		id = MPI::COMM_WORLD.Get_rank ();
	}
	
	messenger::~messenger () {
		MPI::Finalize ();
	}
	
	double& messenger::operator[] (int i) {
		return buffer [i]; 
	}
	
	void messenger::send (double* data, int process, int tag, double weight, int size, int inc) {
		if (weight != 1.0 || inc != 1) {
			if (size > (int) buffer.size ()) {
				buffer.resize (size);
			}
			utils::scale (size, 0.0, &buffer [0]);
			utils::add_scaled (size, weight, data, &buffer [0], inc);
			MPI::COMM_WORLD.Send (&buffer [0], size, MPI::DOUBLE, process, tag);
		} else {
			MPI::COMM_WORLD.Send (data, size, MPI::DOUBLE, process, tag);
		}
	}

	void messenger::recv (double* data, int process, int tag, double weight, int size, int inc) {
		if (weight != 0.0 || inc != 1) {
			if (size > (int) buffer.size ()) {
				buffer.resize (size);
			}
			MPI::COMM_WORLD.Recv (&buffer [0], size, MPI::DOUBLE, process, tag);
			utils::scale (size, weight, data, inc);
			utils::add_scaled (size, 1.0, &buffer [0], data, 1, inc);
		} else {
			MPI::COMM_WORLD.Recv (data, size, MPI::DOUBLE, process, tag);
		}
	}

	void messenger::send (int* data, int process, int tag, int weight, int size, int inc) {
		TRACE ("Sending...");
		if (weight != 1 || inc != 1) {
			if (size > (int) int_buffer.size ()) {
				int_buffer.resize (size);
			}
			for (int i = 0; i < size; ++i) {
				int_buffer [i] = weight * data [i];
			}
			MPI::COMM_WORLD.Send (&int_buffer [0], size, MPI::INT, process, tag);
		} else {
			MPI::COMM_WORLD.Send (data, size, MPI::INT, process, tag);
		}
	}

	void messenger::recv (int* data, int process, int tag, int weight, int size, int inc) {
		TRACE ("Recving...");
		if (weight != 0 || inc != 1) {
			if (size > (int) int_buffer.size ()) {
				int_buffer.resize (size);
			}
			MPI::COMM_WORLD.Recv (&int_buffer [0], size, MPI::INT, process, tag);
			for (int i = 0; i < size; ++i) {
				data [i] = weight * data [i] + int_buffer [i];
			}
		} else {
			MPI::COMM_WORLD.Recv (data, size, MPI::INT, process, tag);
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
} /* utils */