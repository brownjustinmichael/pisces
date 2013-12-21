///////////////////////////////////////////////////////////////////////////////////////
//
//		! \file advection_one_d.h
//		File type: header
//		Author: Ryan Moll (ryan.d.moll@gmail.com)
//		Date created: April 4, 2013
//		Description: Contains methods for explicit handling of spatial derivatives
//
///////////////////////////////////////////////////////////////////////////////////////

#ifndef advection_H
#define advection_H

#include <memory>
#include <vector>
#include "../bases/grid.hpp"
#include "plan_one_d.hpp"

namespace one_d
{
	template <class datatype>
	class advec : public explicit_plan <datatype>
	{
	private:
		using explicit_plan <datatype>::n;
		using explicit_plan <datatype>::data_in;
		using explicit_plan <datatype>::data_out;
	
		datatype c;
		std::vector<datatype> fac;

	public:

		advec (bases::solver <datatype> &i_solver, datatype i_c);	//constuctor initializes private members to point to input and output vectors
		
		virtual ~advec () {
			// printf ("Destroying one_d advection\n");
		}
		
		void execute (int &element_flags);
	};
}

#endif