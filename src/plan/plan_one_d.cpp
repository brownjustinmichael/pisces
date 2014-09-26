/*!**********************************************************************
 * \file plan_one_d.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-11.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "plan_one_d.hpp"

namespace one_d
{
	template class explicit_plan <double>;
	template class explicit_plan <float>;
	
	template class implicit_plan <double>;
	template class implicit_plan <float>;
} /* one_d */