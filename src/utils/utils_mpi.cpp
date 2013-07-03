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

namespace utils
{
	messenger::messenger (int* argc, char*** argv) {
		MPI::Init (*argc, *argv);
		np = MPI::COMM_WORLD.Get_size ();
		id = MPI::COMM_WORLD.Get_rank ();
	}
	
	messenger::~messenger () {
		MTRACE ("Calling destructor");
		MPI::Finalize ();
	}
	
	double& messenger::operator[] (int i) {
		return buffer [i]; 
	}
	
	void messenger::send (double* data, int process, int tag, int size) {
		MPI::COMM_WORLD.Send (data, size, MPI::DOUBLE, process, tag);
	}

	void messenger::recv (int process, int tag, int size) {
		if (size > (int) buffer.size ()) {
			buffer.resize (size);
		}
		MPI::COMM_WORLD.Recv (&buffer [0], size, MPI::DOUBLE, process, tag);
	}
	
	void messenger::min (double* data) {
		if (np != 1) {
			if (np > (int) buffer.size ()) {
				buffer.resize (np);
			}
			MPI::COMM_WORLD.Gather (data, 1, MPI_DOUBLE, &buffer [0], 1, MPI_DOUBLE, 0);
			if (id == 0) {
				*data = *std::min_element (buffer.begin (), buffer.end ());
			}
			MPI::COMM_WORLD.Bcast (data, 1, MPI_DOUBLE, 0);
		}
	}
} /* utils */