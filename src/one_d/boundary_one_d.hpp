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
		boundary (element* i_element_plus, bool i_plus_n, double i_alpha_plus = 0.0, element* i_element_minus = NULL, bool i_minus_n, double i_alpha_minus = 0.0, int *i_flags_ptr = NULL, int i_logger = -1) : bases::boundary (i_flags_ptr, i_logger) {
			alpha_plus = i_alpha_plus;
			alpha_minus = i_alpha_minus;
			element_plus = i_element_plus;
			element_minus = i_element_minus;
			
			plus_n = i_plus_n;
			minus_n = i_minus_n;
			
			if (plus_n) {
				element_plus->flags |= fixed_n;
				index_plus = (*element_plus).n - 1;
			} else {
				element_plus->flags |= fixed_0;
				index_plus = 0;
			}
			
			if (minus_n) {
				element_minus->flags |= fixed_n;
				index_minus = (*element_minus).n - 1;
			} else {
				element_minus->flags |= fixed_0;
				index_minus = 0;
			}
		}
	
		virtual ~boundary () {}
	
		/*!*******************************************************************
		 * \copydoc bases::boundary::execute ()
		 *********************************************************************/
		inline virtual void execute () {
			TRACE (logger, "Executing...");
		
			for (element_plus->iterator iter = element_plus->begin (); iter != element_plus->end (); ++iter) {
				if (!element_minus) {
					(*element_plus) (iter->first, index_plus) = alpha_plus;
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
