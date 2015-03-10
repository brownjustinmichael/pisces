/*!**********************************************************************
 * \file element_two_d.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-12.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "implemented_element.hpp"

namespace pisces
{
	// template <class datatype>
	// int implemented_element <datatype>::mode = mode_flag;
	//
	// template class element <double>;
		
	template <class datatype>
	int implemented_element <datatype>::mode = grids::mode_flag;

	template class implemented_element <double>;
} /* pisces */