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
#include "advection_one_d.hpp"
#include "plan_one_d.hpp"
#include "../utils/utils.hpp"
#include "../bases/element.hpp"

namespace one_d
{
	template <class datatype>
	advec <datatype>::advec (bases::grid <datatype> &i_grid, datatype i_c, datatype* i_data_in, datatype* i_data_out) : 
	explicit_plan <datatype> (i_grid, i_data_in, i_data_out)
	{
		TRACE ("Instantiating...");
		datatype pi = std::acos(-1.0);
		c = i_c;
		
		fac.resize(n,0.0);

		fac [0] = std::sin(((n - 1.5)*pi)/(n-1));
		for (int i = 1; i < n; i++)
		{
			fac [i] = std::sin(((n - i - 1)*pi)/(n-1));
		}
		fac [n - 1] = std::sin((.5*pi)/(n-1));
		TRACE ("Instantiated.");
	}

	template <class datatype>
	void advec <datatype>::execute(int &element_flags)
	{
		datatype scalar = -c;
		
		// if (!(*flags_ptr | transformed)) {
		// 		FATAL ("Linear Advection attempted in physical space.")
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
		
		data_out [0] += (scalar/fac [0])*(data_in [1] - data_in [0])*data_in [0];
		for (int i = 1; i < (n - 1); i++)
		{
			data_out [i] += .5*(scalar/fac [i])*(data_in [i + 1] - data_in [i - 1])*data_in [i];
		}
		data_out [n - 1] += (scalar/fac [n - 1])*(data_in [n - 1] - data_in [n - 2])*data_in [n - 1];
	}
	
	template class advec <double>;
	template class advec <float>;
}
