///////////////////////////////////////////////////////////////////////////////////////
//
//		! \file advection_one_d.cpp
//		File type: source
//		Author: Ryan Moll (ryan.d.moll@gmail.com)
//		Date created: April 4, 2013
//		Description: Contains methods for explicit handling of spatial derivatives
//
///////////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "../bases/boundary.hpp"
#include "advection_one_d.h"
#include "../bases/plan.hpp"
#include "../utils/utils.hpp"

namespace one_d
{
	
	advec::advec (int i_n, double i_c, int i_name_in, int i_name_out, std::shared_ptr<bases::collocation_grid> i_grid) : bases::explicit_plan (i_n, i_name_in, i_name_out)
	{
		MTRACE ("Instantiating...");
		grid = i_grid;
		c = i_c;
		MTRACE ("Instantiated.");
	}

	void advec::execute()
	{
		bases::plan::execute ();

		if (*flags_ptr & linked_0) {
			data_out [0] -= c * utils::dot (n, grid->get_data (1), data_in, n);
		}

		utils::matrix_vector_multiply (n - 2, n, -c, grid->get_data (1) + 1, data_in, 1.0, data_out + 1, n);
		
		if (*flags_ptr & linked_n) {
			data_out [n - 1] -= c * utils::dot (n, grid->get_data (1) + n - 1, data_in, n);
		}
	}
}
