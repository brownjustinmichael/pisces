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
	template <class datatype>
	std::shared_ptr <typename plan <datatype>::factory> advec (grids::variable &vel_n, grids::variable &vel_m) {
		return std::shared_ptr <typename real_plan <datatype>::factory> (new typename advection::uniform <datatype>::factory (vel_n, vel_m));
	}
} /* plans */

#endif /* end of include guard: ADVECTION_HPP_5C2EB309 */
