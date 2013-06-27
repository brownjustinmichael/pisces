/*!***********************************************************************
 * \file boundary_one_d.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-05-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "boundary_one_d.hpp"
#include "mpi.h"

namespace one_d
{
	void mpi_boundary::send (int name) {
		TRACE (logger, "Sending Tag: " << ext_send << " to Process: " << process);
		buffer = (&((*element_ptr) [name])) [index];
		MPI::COMM_WORLD.Send (&buffer, 1, MPI::DOUBLE, process, ext_send);
		TRACE (logger, "Sent.")
	}

	void mpi_boundary::recv (int name) {
		TRACE (logger, "Recving Tag: " << ext_recv << " from Process: " << process);
		MPI::COMM_WORLD.Recv (&buffer, 1, MPI::DOUBLE, process, ext_recv);
		(&((*element_ptr) [name])) [index] = 0.5 * (&((*element_ptr) [name])) [index] + 0.5 * buffer;
		TRACE (logger, "Recved.")
	}
} /* one_d */
