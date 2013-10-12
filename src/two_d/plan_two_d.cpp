/*!**********************************************************************
 * \file plan_two_d.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-11.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "plan_two_d.hpp"

namespace two_d
{
	template <class datatype>
	explicit_plan <datatype>::explicit_plan (int i_n, int i_m, datatype* i_data_in, datatype* i_data_out) : 
	n (i_n),
	m (i_m),
	data_in (i_data_in) {
		TRACE ("Instantiating...");
		if (i_data_out == NULL) {
			data_out = data_in;
		} else {
			data_out = i_data_out;
		}
		TRACE ("Instantiated.");
	}
	
	template class explicit_plan <double>;
	template class explicit_plan <float>;
	
	template class implicit_plan <double>;
	template class implicit_plan <float>;
} /* two_d */