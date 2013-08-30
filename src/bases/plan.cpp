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
	explicit_plan <datatype>::explicit_plan (int i_n, datatype* i_data_in, datatype* i_data_out) : 
	n (i_n),
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
} /* bases */