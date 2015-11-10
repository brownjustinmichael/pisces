/*!**********************************************************************
 * \file boundary.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-07-14.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "boundary.hpp"

namespace boundaries
{
	void boundary_match (int m, std::shared_ptr <boundaries::boundary> &boundary_0, std::shared_ptr <boundaries::boundary> &boundary_n, double *data, int excess_0, int excess_n) {
		if (boundary_n) {
			boundary_n->send (data + m - 2 * excess_n - 1, m, excess_n);
		}
		if (boundary_0) {
			boundary_0->receive (data, m, excess_0, 0.0);
			boundary_0->send (data + excess_0, m, excess_0 + 1);
		}
		if (boundary_n) {
			boundary_n->receive (data + m - excess_n - 1, m, excess_n + 1, 0.0);
		}
	}
} /* boundaries */
