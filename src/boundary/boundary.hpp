/*!***********************************************************************
 * \file boundary.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BOUNDARY_HPP_CP5DR4CP
#define BOUNDARY_HPP_CP5DR4CP

#include "../plan.hpp"

namespace boundary
{
	enum flags {
		fixed_upper = 0x01,
		fixed_lower = 0x02
	};
} /* boundary */

#endif /* end of include guard: BOUNDARY_HPP_CP5DR4CP */
