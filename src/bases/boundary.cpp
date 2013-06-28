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
	boundary::boundary (element* i_element_ptr, int i_edge, int i_process) : plan (i_element_ptr) {
		edge = i_edge;
		process = i_process;
		element_ptr->get_boundary_info (edge, index, increment);
		*flags_ptr |= edge;
	}
} /* bases */
