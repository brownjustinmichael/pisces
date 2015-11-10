/*!**********************************************************************
 * \file plans/advection.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "advection.hpp"

namespace plans
{
	std::shared_ptr <plan::factory> advec (grids::variable &vel_n, grids::variable &vel_m) {
		return std::shared_ptr <real_plan::factory> (new advection::uniform::factory (vel_n, vel_m));
	}
} /* plans */

