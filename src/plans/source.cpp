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
	std::shared_ptr <typename plan <double>::factory> src (grids::variable &data_source, bool dealias) {
		if (data_source.get_states () > 1) {
			return std::shared_ptr <typename explicit_plan <double>::factory> (new typename source::uniform <double>::factory (data_source, dealias, 1.0));
		} else {
			return std::shared_ptr <typename real_plan <double>::factory> (new typename source::uniform_real <double>::factory (data_source, 1.0));
		}
	}

	std::shared_ptr <typename plan <double>::factory> constant (double coeff) {
		return std::shared_ptr <typename explicit_plan <double>::factory> (new typename source::constant <double>::factory (coeff));
	}
} /* plans */
