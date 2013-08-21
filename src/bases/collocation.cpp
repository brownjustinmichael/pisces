/*!***********************************************************************
 * \file bases/collocation.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef COLLOCATION_CPP_HV4P0UOP
#define COLLOCATION_CPP_HV4P0UOP

#include "../config.hpp"
#include "collocation.hpp"

namespace bases
{
	template class collocation_grid <double>;
	template class collocation_grid <float>;
} /* bases */

#endif /* end of include guard: COLLOCATION_CPP_HV4P0UOP */
