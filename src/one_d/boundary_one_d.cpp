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
	void link_boundary::send () {
		int j = 0;
		TRACE (logger, "Sending Tag: " << ext_send << " to Process: " << process);
		for (bases::element::iterator iter = (*element_ptr).begin (); iter != (*element_ptr).end (); ++iter) {
			for (int i = 0; i < n; ++i) {
				to_send [i + n * j] = send_buffer [*iter] [i] = (&((*element_ptr) [*iter])) [index + i * increment];
			}
			++j;
		}
		MPI::COMM_WORLD.Send (&to_send [0], j * n, MPI::DOUBLE, process, ext_send);
		TRACE (logger, "Sent.")
	}

	void link_boundary::recv () {
		int j = 0;
		TRACE (logger, "Recving Tag: " << ext_recv << " from Process: " << process << " ptr " << &to_recv [0] << " " << &to_recv [n * j]);
		for (bases::element::iterator iter = (*element_ptr).begin (); iter != (*element_ptr).end (); ++iter) {
			++j;
		}
		TRACE (logger, "here")
		MPI::COMM_WORLD.Recv (&to_recv [0], j * n, MPI::DOUBLE, process, ext_recv);
		j = 0;
		TRACE (logger, "here")
		for (bases::element::iterator iter = (*element_ptr).begin (); iter != (*element_ptr).end (); ++iter) {
			for (int i = 0; i < n; ++i) {
				TRACE (logger, "ptr " << &recv_buffer [*iter] [0]);
				recv_buffer [*iter] [i] = to_recv [i + n * j];
			}
			++j;
		}
		TRACE (logger, "Recved.")
	}
} /* one_d */
