/*!**********************************************************************
 * \file real_plan.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "real_plan.hpp"

namespace plans
{
	std::shared_ptr <typename compound_plan <double>::factory> operator* (std::shared_ptr <typename plan <double>::factory> i_factory1, std::shared_ptr <typename plan <double>::factory> i_factory2) {
		std::shared_ptr <typename compound_plan <double>::factory> plan = typename std::shared_ptr <typename compound_plan <double>::factory> (new typename compound_plan <double>::factory ());

		plan->add_plan (i_factory1, compound_plan <double>::mult);
		plan->add_plan (i_factory2, compound_plan <double>::mult);
		return plan;
	}

	std::shared_ptr <typename compound_plan <double>::factory> operator/ (std::shared_ptr <typename plan <double>::factory> i_factory1, std::shared_ptr <typename plan <double>::factory> i_factory2) {
		std::shared_ptr <typename compound_plan <double>::factory> plan = std::shared_ptr <typename compound_plan <double>::factory> (new typename compound_plan <double>::factory ());

		plan->add_plan (i_factory1, compound_plan <double>::mult);
		plan->add_plan (i_factory2, compound_plan <double>::div);
		return plan;
	}
} /* plans */
