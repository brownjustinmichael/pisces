/*!**********************************************************************
 * \file data.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-10-09.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DATA_CPP_A42CF84C
#define DATA_CPP_A42CF84C

#include "data.hpp"

namespace data
{
	template <class datatype>
	int data <datatype>::mode = plans::mode_flag;
	
	template class data <double>;
} /* data */

#endif /* end of include guard: DATA_CPP_A42CF84C */
