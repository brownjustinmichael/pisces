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
	template <class datatype>
	std::shared_ptr <typename plan <datatype>::factory> advec (grids::variable <datatype> &vel_n, grids::variable <datatype> &vel_m, datatype coeff = 1.0) {
		return std::shared_ptr <typename real_plan <datatype>::factory> (new typename advection::uniform <datatype>::factory (vel_n, vel_m, coeff));
	}
} /* plans */

#endif /* end of include guard: ADVECTION_HPP_5C2EB309 */
