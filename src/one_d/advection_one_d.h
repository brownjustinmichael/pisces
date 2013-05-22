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
	class advec : public bases::explicit_plan
	{
	private:
	
		double c;
		std::shared_ptr<bases::collocation_grid> grid;

	public:

		advec (int i_n, double i_c, int i_name_in, int i_name_out, std::shared_ptr<bases::collocation_grid> i_grid);	//constuctor initializes private members to point to input and output vectors
		
		virtual ~advec () {}
		
		void execute ();
	};
}

#endif