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
		if (ext_index) {
			ext_index = ext_element_ptr->get_boundary_index (ext_edge);
		} else {
			*flags_ptr |= edge;
		}
		MTRACE ("Associated." << (*element_ptr) [0]);
	}

	void boundary::execute () {
		plan::execute ();
	
		for (int i = 0; i < element_ptr->n_explicit_grid_plans; ++i) {
			element_ptr->explicit_grid_plans [i]->boundary (index, ext_element_ptr, ext_index);
		}

		for (int i = 0; i < element_ptr->n_explicit_space_plans; ++i) {
			element_ptr->explicit_space_plans [i]->boundary (index, ext_element_ptr, ext_index);
		}

		if (!(*flags_ptr & factorized)) {
			for (int i = 0; i < element_ptr->n_implicit_plans; ++i) {
				element_ptr->implicit_plans [i]->boundary (index, ext_element_ptr, ext_index);
			}
		}
	}
} /* bases */
