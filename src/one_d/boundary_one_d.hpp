/*!***********************************************************************
 * \file one_d/boundary.hpp
 * Spectral Element
 * 
 * This file provides the implementation for the boundary base class. 
 * New realizations of boundary conditions should include this file to 
 * inherit from the boundary class.
 * 
 * Created by Justin Brown on 2013-04-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BOUNDARY_HPP_CP5DR4CP
#define BOUNDARY_HPP_CP5DR4CP

#include <memory>
#include "../bases/boundary.hpp"
#include "element.hpp"
#include "../config.hpp"

namespace one_d
{
	/*!*******************************************************************
	 * \brief An implementation of the boundary class in 1D
	 * 
	 * This class handles the single point boundary by taking a weighted 
	 * average of the edges of two elements. It can also be used to zero 
	 * a boundary.
	 *********************************************************************/
	class boundary : public bases::boundary
	{
	public:
		/*!*******************************************************************
		 * \copydoc bases::boundary::boundary ()
		 *********************************************************************/
		boundary (int i_edge = fixed_0, double i_alpha = 0.0, bases::element* i_ext_element_ptr = NULL, int i_ext_edge = fixed_0, double i_ext_alpha = 0.0) : bases::boundary (i_edge, i_alpha, i_ext_element_ptr, i_ext_edge, i_ext_alpha) {}
	
		virtual ~boundary () {}
	
		/*!*******************************************************************
		 * \copydoc bases::boundary::execute ()
		 *********************************************************************/
		inline virtual void execute () {
			TRACE (logger, "Executing...");
		
			for (element_plus->iterator iter = element_plus->begin (); iter != element_plus->end (); ++iter) {
				if (!ext_element_ptr) {
					(*ext_element_ptr) (iter->first, index_plus) = alpha_plus;
				} else if (!element_plus) {
					(*element_minus) (iter->first, index_minus) = alpha_minus;
				} else {
					*data_minus = alpha_plus * *data_plus + alpha_minus * *data_minus;
					if (!one_way) {
						*data_plus = *data_minus;					
					}
				}
			}
			
			/*
				TODO Make implementation more general
			*/
			

		
			TRACE (logger, "executed.")
		}
	private:
		double alpha_plus; //!< A double coefficient for the contribution from the positive boundary
		double alpha_minus; //!< A double coefficient for the contribution from the negative boudary
		element* element_plus; //!< A pointer to the double first element of the positive boundary
		element* element_minus; //!< A pointer to the double first element of the positive boundary
		bool plus_n;
		bool minus_n;
		int index_plus;
		int index_minus;
	};
} /* one_d */


#endif /* end of include guard: BOUNDARY_HPP_CP5DR4CP */
