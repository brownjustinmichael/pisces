/*!***********************************************************************
 * \file plan.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-30.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "../config.hpp"
#include "plan.hpp"
#include "element.hpp"

namespace bases
{
	void plan::associate (element* i_element_ptr) {
		MTRACE ("Associating...");
		element_ptr = i_element_ptr;
		logger = element_ptr->logger;
		flags_ptr = &(element_ptr->flags);
		MTRACE ("Associated.")
	}
	
	void explicit_plan::associate (element* i_element_ptr) {
		plan::associate (i_element_ptr);
		data_in = &((*i_element_ptr) [name_in]);
		data_out = &((*i_element_ptr) [name_out]);
	}
} /* bases */