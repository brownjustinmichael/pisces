/*!***********************************************************************
 * \file boundary.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-05-01.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "plan.hpp"
#include "boundary.hpp"
#include "element.hpp"

namespace bases
{
	void boundary::associate (bases::element* i_element_ptr) {
		bases::plan::associate (i_element_ptr);
		MTRACE ("Associating...");
		index = element_ptr->get_boundary_index (edge);
		if (ext_element_ptr) {
			ext_index = ext_element_ptr->get_boundary_index (ext_edge);
		} else {
			*flags_ptr |= edge;
		}
		MTRACE ("Associated." << (*element_ptr) [0]);
	}

	void boundary::send () {
		int j;
		j = 0;
		send_buffer.resize (n * (*element_ptr).names.size ());
		for (bases::element::iterator iter = (*element_ptr).begin (); iter != (*element_ptr).end (); ++iter) {
			for (int i = 0; i < n; ++i) {
				send_buffer [i + n * j] = (&((*element_ptr) [*iter])) [index + i * to_next];
			}
			++j;
		}
	}
	
	void boundary::recv () {
		int j;
		j = 0;
		recv_buffer.resize (n * (*element_ptr).names.size ());
		for (bases::element::iterator iter = (*element_ptr).begin (); iter != (*element_ptr).end (); ++iter) {
			for (int i = 0; i < n; ++i) {
				recv_buffer [i + n * j] = ext_buffer [i + n * j];
			}
			++j;
		}
	}
	
	void active_boundary::execute () {
		plan::execute ();
		
		for (int i = 0; i < n_plans; ++i) {
			plans [i]->boundary (index, ext_element_ptr, ext_index);
		}
	}
} /* bases */
