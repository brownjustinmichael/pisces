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
	plan::factory_container diff (double alpha) {
		return implicit_plan::factory_container (std::shared_ptr <implicit_plan::factory> (new diffusion::vertical::factory (1.0, alpha))) + 
		implicit_plan::factory_container (std::shared_ptr <implicit_plan::factory> (new diffusion::horizontal::factory (1.0, alpha)));
	}

	plan::factory_container bg_diff (double *i_diffusion, double alpha) {
		return implicit_plan::factory_container (std::shared_ptr <plan::factory> (new diffusion::background_vertical::factory (alpha, i_diffusion))) + 
		implicit_plan::factory_container (std::shared_ptr <plan::factory> (new diffusion::background_horizontal::factory (alpha, i_diffusion)));
	}

	plan::factory_container density_diff (grids::variable &density, double alpha) {
		return diff (alpha) +
		plan::factory_container (std::shared_ptr <plan::factory> (new diffusion::variable_diffusion::factory (density, 1.0)));
	}

	std::shared_ptr <plan::factory> horizontal_stress (grids::variable &density, grids::variable &data_other) {
		return std::shared_ptr <explicit_plan::factory> (new diffusion::horizontal_stress::factory (density, data_other, 1.0));
	}

	std::shared_ptr <plan::factory> vertical_stress (grids::variable &density, grids::variable &data_other) {
		return std::shared_ptr <explicit_plan::factory> (new diffusion::vertical_stress::factory (density, data_other, 1.0));
	}
} /* plans */
