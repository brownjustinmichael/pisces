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
	template <class datatype>
	plan <datatype>::plan (element <datatype>* i_element_ptr, int i_flags) :
	element_ptr (i_element_ptr),
	flags (i_flags) {
		TRACE ("Instantiating...");
		flags_ptr = &(element_ptr->flags);
		messenger_ptr = element_ptr->messenger_ptr;
		TRACE ("Instantiated.");
	}
	
	template <class datatype>
	explicit_plan <datatype>::explicit_plan (element <datatype>* i_element_ptr, int i_n, int i_name_in, int i_name_out, int i_flags) : 
	plan <datatype> (i_element_ptr, i_flags),
	n (i_n) {
		TRACE ("Instantiating...");
		data_in = &((*i_element_ptr) [i_name_in]);
		if (i_name_out == null) {
			data_out = &((*i_element_ptr) [i_name_in]);
		} else {
			data_out = &((*i_element_ptr) [i_name_out]);
		}
		TRACE ("Instantiated.");
	}
	
	template <class datatype>
	void explicit_plan <datatype>::execute () {
		plan <datatype>::execute ();
	}
	
	template <class datatype>
	void implicit_plan <datatype>::execute () {
		plan <datatype>::execute ();
	}
	
	template class plan <double>;
	template class plan <float>;
	
	template class explicit_plan <double>;
	template class explicit_plan <float>;
	
	template class implicit_plan <double>;
	template class implicit_plan <float>;
} /* bases */