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
#include "advection_one_d.h"
#include "../bases/transform.hpp"
#include "../bases/plan.hpp"
#include "../utils/utils.hpp"
#include "../bases/element.hpp"

namespace one_d
{
	advec::advec (bases::element* i_element_ptr, int i_n, double i_c, int i_name_in, int i_name_out, std::shared_ptr<bases::collocation_grid> i_grid) : bases::explicit_plan (i_element_ptr, i_n, i_name_in, i_name_out)
	{
		MTRACE ("Instantiating...");
		grid = i_grid;
		double pi = std::acos(-1.0);
		c = i_c;
		
		fac.resize(n,0.0);

		fac [0] = std::sin(((n - 1.5)*pi)/(n-1));
		for (int i = 1; i < n; i++)
		{
			fac [i] = std::sin(((n - i - 1)*pi)/(n-1));
		}
		fac [n - 1] = std::sin((.5*pi)/(n-1));
		MTRACE ("Instantiated.");
	}

	void advec::execute()
	{
		double scalar = -c;
		
		bases::plan::execute ();

		// if (!(*flags_ptr | transformed)) {
		// 		MFATAL ("Linear Advection attempted in physical space.")
		// 		throw 0;
		// 	}
		//  
		// if (*flags_ptr & linked_0) {
		// 	data_out [0] -= c * utils::dot (n, grid->get_data (1), data_in, n);
		// }
		// 
		// utils::matrix_vector_multiply (n - 2, n, -c, grid->get_data (1) + 1, data_in, 1.0, data_out + 1, n);
		// 
		// if (*flags_ptr & linked_n) {
		// 	data_out [n - 1] -= c * utils::dot (n, grid->get_data (1) + n - 1, data_in, n);
		// }
		
		if (*flags_ptr & transformed) {
			MFATAL ("Nonlinear advection attempted in Chebyshev space.")
			throw 0;
		}
		
		data_out [0] += (scalar/fac [0])*(data_in [1] - data_in [0])*data_in [0];
		for (int i = 1; i < (n - 1); i++)
		{
			data_out [i] += .5*(scalar/fac [i])*(data_in [i + 1] - data_in [i - 1])*data_in [i];
		}
		data_out [n - 1] += (scalar/fac [n - 1])*(data_in [n - 1] - data_in [n - 2])*data_in [n - 1];
	}
}
