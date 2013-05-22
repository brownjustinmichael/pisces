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
	void boundary::associate (element* i_element_ptr) {
		plan::associate (i_element_ptr);
		MTRACE ("Associating...");
		element_ptr->get_boundary_info (edge, index, increment);
		*flags_ptr |= edge;
		MTRACE ("Associated.");
	}
} /* bases */
