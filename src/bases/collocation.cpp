/*!***********************************************************************
 * \file bases/collocation.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef COLLOCATION_CPP_HV4P0UOP
#define COLLOCATION_CPP_HV4P0UOP

#include "../config.hpp"
#include "collocation.hpp"

namespace bases
{
	collocation_grid::collocation_grid (int i_derivs, int i_rows, int i_cols, int i_logger) {
		logger = i_logger;
		rows = i_rows;
		cols = i_cols;
		derivs = i_derivs;

		TRACE (logger, "Instantiating...")

		data.resize (derivs);

		for (int i = 0; i < i_derivs; ++i) {
			data [i].resize (i_rows * i_cols);
		}

		TRACE (logger, "Instantiated...")
	}
} /* bases */

#endif /* end of include guard: COLLOCATION_CPP_HV4P0UOP */
