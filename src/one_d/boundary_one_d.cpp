/*!***********************************************************************
 * \file boundary_one_d.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-05-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "boundary_one_d.hpp"

namespace one_d
{
	void link_boundary::send () {
		int j = 0;
		for (bases::element::iterator iter = (*element_ptr).begin (); iter != (*element_ptr).end (); ++iter) {
			for (int i = 0; i < n; ++i) {
				ext_send [i + n * j] = send_buffer [*iter] [i] = (&((*element_ptr) [*iter])) [index + i * increment];
				MDEBUG ("" << (i + n * j) << " " << ext_send [i+n*j]);
			}
			++j;
		}
	}

	void link_boundary::recv () {
		int j = 0;
		for (bases::element::iterator iter = (*element_ptr).begin (); iter != (*element_ptr).end (); ++iter) {
			for (int i = 0; i < n; ++i) {
				recv_buffer [*iter] [i] = ext_recv [i + n * j];
				MDEBUG ("" << (i + n * j) << " " << ext_recv [i+n*j]);
			}
			++j;
		}
	}
} /* one_d */
