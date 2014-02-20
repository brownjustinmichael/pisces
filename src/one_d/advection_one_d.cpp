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

namespace one_d
{
	template <class datatype>
	advection <datatype>::advection (bases::grid <datatype> &i_grid, datatype i_coeff, datatype *i_data_in, datatype *i_data_out, int *i_element_flags, int *i_component_flags) :
	explicit_plan <datatype> (i_grid, i_data_in, i_data_out, i_element_flags, i_component_flags),
	coeff (-i_coeff),
	position_ptr (&i_grid [0]) {}
	
	template <class datatype>
	advection <datatype>::advection (bases::solver <datatype> &i_solver, datatype i_coeff) : 
	explicit_plan <datatype> (i_solver),
	coeff (-i_coeff),
	position_ptr (&grid [0]) {}

	template <class datatype>
	void advection <datatype>::execute()
	{
		data_out [0] += coeff * (data_in [1] - data_in [0]) / (position_ptr [1] - position_ptr [0]) * data_in [0];
		for (int i = 1; i < n - 1; ++i)
		{
			data_out [i] += coeff * (data_in [i + 1] - data_in [i - 1]) / (position_ptr [i + 1] - position_ptr [i - 1]) * data_in [i];
		}
		data_out [n - 1] += coeff * (data_in [n - 1] - data_in [n - 2]) / (position_ptr [n - 1] - position_ptr [n - 2]) * data_in [n - 1];
	}
	
	template class advection <double>;
	template class advection <float>;
}
