/*!**********************************************************************
 * \file plans/advection.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ADVECTION_HPP_5C2EB309
#define ADVECTION_HPP_5C2EB309

#include "plans-advection/advection.hpp"

#include "real_plan.hpp"

namespace plans
{
	/**
	 * @brief Shorthand to generate an advection plan
	 * 
	 * @param vel_n A reference to the x component of the velocity
	 * @param vel_m A reference to the z component of the velocity
	 * 
	 * @return [description]
	 */
	std::shared_ptr <plan::factory> advec (grids::variable &vel_n, grids::variable &vel_m);
} /* plans */

#endif /* end of include guard: ADVECTION_HPP_5C2EB309 */
