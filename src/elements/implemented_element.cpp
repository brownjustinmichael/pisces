/*!**********************************************************************
 * \file implemented_element.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-12.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "implemented_element.hpp"

namespace pisces
{
	// template <class datatype>
	// int implemented_element::mode = mode_flag;
	//
	// template class element;

	std::map <std::string, implemented_element::element_function> & implemented_element::registry()
	{
	    static std::map <std::string, implemented_element::element_function> impl;
	    return impl;
	}

	int implemented_element::mode = grids::mode_flag;
} /* pisces */