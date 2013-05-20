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
		for (bases::element::iterator iter = (*element_ptr).begin (); iter != (*element_ptr).end (); ++iter) {
			for (int i = 0; i < n; ++i) {
				to_send [i + n * j] = send_buffer [*iter] [i] = (&((*element_ptr) [*iter])) [index + i * increment];
				MDEBUG ("send " << *iter << " " << i << " " << send_buffer [*iter] [i]);
			}
			++j;
		}
		MPI::COMM_WORLD.Send (&to_send [0], (j + 1) * n, MPI::REAL, 0, ext_send);
	}

	void link_boundary::recv () {
		int j = 0;
		for (bases::element::iterator iter = (*element_ptr).begin (); iter != (*element_ptr).end (); ++iter) {
			++j;
		}
		MPI::COMM_WORLD.Recv (&to_recv [0], (j + 1) * n, MPI::REAL, 0, ext_recv);
		j = 0;
		for (bases::element::iterator iter = (*element_ptr).begin (); iter != (*element_ptr).end (); ++iter) {
			for (int i = 0; i < n; ++i) {
				recv_buffer [*iter] [i] = to_recv [i + n * j];
				MDEBUG ("recv " << *iter << " " << i << " " << recv_buffer [*iter] [i]);
			}
			++j;
		}
	}
} /* one_d */
