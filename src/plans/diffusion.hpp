/*!**********************************************************************
 * \file diffusion.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DIFFUSION_HPP_C1935DDC
#define DIFFUSION_HPP_C1935DDC

#include "plans-diffusion/horizontal_diffusion.hpp"
#include "plans-diffusion/vertical_diffusion.hpp"
#include "plans-diffusion/variable_diffusion.hpp"

#include "implicit_plan.hpp"
#include "io/parameters.hpp"

namespace plans
{
	template <class datatype>
	typename implicit_plan <datatype>::factory_container diff (datatype coeff = 1.0) {
		return typename implicit_plan <datatype>::factory_container (std::shared_ptr <typename implicit_plan <datatype>::factory> (new typename diffusion::vertical <datatype>::factory (coeff))) + typename implicit_plan <datatype>::factory_container (std::shared_ptr <typename implicit_plan <datatype>::factory> (new typename diffusion::horizontal <datatype>::factory (coeff)));
	}
} /* plans */

#endif /* end of include guard: DIFFUSION_HPP_C1935DDC */
