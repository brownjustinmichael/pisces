/*!***********************************************************************
 * \file bases/message.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "messenger.hpp"
#include "../config.hpp"
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
	
	void messenger::send (double* data, int process, int tag, int size) {
		MPI::COMM_WORLD.Send (data, size, MPI::DOUBLE, process, tag);
	}

	void messenger::recv (double* data, int process, int tag, int size) {
		MPI::COMM_WORLD.Recv (data, size, MPI::DOUBLE, process, tag);
	}

} /* utils */