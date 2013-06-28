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
		MTRACE ("Instantiating...");
		element_ptr = i_element_ptr;
		logger = element_ptr->logger;
		flags_ptr = &(element_ptr->flags);
		MTRACE ("Instantiated.");
	}
	
	explicit_plan::explicit_plan (element* i_element_ptr, int i_n, int i_name_in, int i_name_out) : plan (i_element_ptr) {
		MTRACE ("Instantiating...");
		n = i_n;
		data_in = &((*i_element_ptr) [i_name_in]);
		if (i_name_out == null) {
			data_out = &((*i_element_ptr) [i_name_in]);
		} else {
			data_out = &((*i_element_ptr) [i_name_out]);
		}
		MTRACE ("Instantiated.");
	}
	
	void explicit_plan::execute () {
		plan::execute ();
		TRACE (logger, "Executing...");
		if (!(*flags_ptr & explicit_started)) {
			element_ptr->explicit_reset ();
			*flags_ptr |= explicit_started;
		}
		TRACE (logger, "Executed.");
	}
	
	void implicit_plan::execute () {
		plan::execute ();
		TRACE (logger, "Executing...");
		if (!(*flags_ptr & implicit_started)) {
			element_ptr->implicit_reset ();
			*flags_ptr |= implicit_started;
		}
		TRACE (logger, "Executed.");
	}
} /* bases */