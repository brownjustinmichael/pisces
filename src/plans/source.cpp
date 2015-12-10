/*!**********************************************************************
 * \file plans/source.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "plans/source.hpp"

namespace plans
{
	std::shared_ptr <plan::factory> src (grids::variable &data_source, bool dealias) {
		if (data_source.get_states () > 1) {
			return std::shared_ptr <explicit_plan::factory> (new source::uniform::factory (data_source, dealias, 1.0));
		} else {
			return std::shared_ptr <real_plan::factory> (new source::uniform_real::factory (data_source, 1.0));
		}
	}

	std::shared_ptr <plan::factory> grad_x (grids::variable &data_source) {
		return std::shared_ptr <explicit_plan::factory> (new source::uniform_grad_x::factory (data_source, 1.0));
	}

	std::shared_ptr <plan::factory> grad_z (grids::variable &data_source) {
		return std::shared_ptr <explicit_plan::factory> (new source::uniform_grad_z::factory (data_source, 1.0));
	}

	plan::factory_container heat (grids::variable &data_source, grids::variable &data_x, grids::variable &data_z) {
		return plan::factory_container (std::shared_ptr <real_plan::factory> (new source::viscous_heat::factory (data_source, data_x, data_z, 1.0)));
	}

	plan::factory_container diverge (grids::variable &data_source, grids::variable &data_x, grids::variable &data_z) {
		return plan::factory_container (std::shared_ptr <real_plan::factory> (new source::divergence::factory (data_source, data_x, data_z, 1.0)));
	}

	std::shared_ptr <plan::factory> constant (double coeff) {
		return std::shared_ptr <explicit_plan::factory> (new source::constant::factory (coeff));
	}
} /* plans */
