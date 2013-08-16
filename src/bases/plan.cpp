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
	plan::plan (element* i_element_ptr) {
		TRACE ("Instantiating...");
		element_ptr = i_element_ptr;
		flags_ptr = &(element_ptr->flags);
		TRACE ("Instantiated.");
	}
	
	explicit_plan::explicit_plan (element* i_element_ptr, int i_n, int i_name_in, int i_name_out) : plan (i_element_ptr) {
		TRACE ("Instantiating...");
		n = i_n;
		data_in = &((*i_element_ptr) [i_name_in]);
		if (i_name_out == null) {
			data_out = &((*i_element_ptr) [i_name_in]);
		} else {
			data_out = &((*i_element_ptr) [i_name_out]);
		}
		TRACE ("Instantiated.");
	}
	
	void explicit_plan::execute () {
		plan::execute ();
	}
	
	void implicit_plan::execute () {
		plan::execute ();
	}
} /* bases */