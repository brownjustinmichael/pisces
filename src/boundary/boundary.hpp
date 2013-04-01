// 
//! \file boundary.hpp
//  boundary
//  
//  Created by Justin Brown on 2013-04-01.
//  Copyright 2013 Justin Brown. All rights reserved.
// 

#ifndef BOUNDARY_HPP_CP5DR4CP
#define BOUNDARY_HPP_CP5DR4CP

#include "../plan.hpp"

namespace boundary
{
	class fixed_cart_1D : public plan
	{
	private:
		double *element; //!< A pointer to the double data element to fix
		double value; //!< The double value at which the element will remained fixed
	public:
		//! \param i_element A pointer to the double data element to fix
		//! \param i_vaule The double value at which the element will remain fixed
		fixed_cart_1D (double *i_element, double i_value) {element = i_element; value = i_value;}
		virtual ~fixed_cart_1D ();
		
		inline void execute (double timestep) {*element = value;}
	};
} /* boundary */

#endif /* end of include guard: BOUNDARY_HPP_CP5DR4CP */
