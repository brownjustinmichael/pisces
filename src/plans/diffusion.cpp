/*!**********************************************************************
 * \file diffusion.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "plans/diffusion.hpp"

namespace plans
{
	std::shared_ptr <typename plan <double>::factory> horizontal_stress (grids::variable &density, grids::variable &data_other) {
		return std::shared_ptr <typename explicit_plan <double>::factory> (new typename diffusion::horizontal_stress <double>::factory (density, data_other, 1.0));
	}

	std::shared_ptr <typename plan <double>::factory> vertical_stress (grids::variable &density, grids::variable &data_other) {
		return std::shared_ptr <typename explicit_plan <double>::factory> (new typename diffusion::vertical_stress <double>::factory (density, data_other, 1.0));
	}
} /* plans */
