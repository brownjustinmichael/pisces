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
#include "../bases/collocation.hpp"
#include "../bases/plan.hpp"

namespace one_d
{
	template <class datatype>
	class advec : public bases::explicit_plan <datatype>
	{
	private:
		using bases::explicit_plan <datatype>::n;
		using bases::explicit_plan <datatype>::data_in;
		using bases::explicit_plan <datatype>::data_out;
		using bases::explicit_plan <datatype>::flags_ptr;
	
		datatype c;
		std::vector<datatype> fac;
		std::shared_ptr<bases::collocation_grid <datatype>> grid;

	public:

		advec (bases::element <datatype>* i_element_ptr, int i_n, datatype i_c, int i_name_in, int i_name_out, std::shared_ptr<bases::collocation_grid <datatype>> i_grid);	//constuctor initializes private members to point to input and output vectors
		
		virtual ~advec () {}
		
		void execute ();
	};
}

#endif