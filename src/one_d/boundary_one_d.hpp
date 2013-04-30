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
		boundary (double *i_data_plus, double i_alpha_plus = 0.0, double *i_data_minus = NULL, double i_alpha_minus = 0.0, bool i_one_way = false, int *i_flags_ptr = NULL, int i_logger = -1) : bases::boundary (i_data_plus, i_alpha_plus, i_data_minus, i_alpha_minus, i_flags_ptr, i_logger) {one_way = i_one_way;}
		
		boundary (double& i_data_plus, double i_alpha_plus, double& i_data_minus, double i_alpha_minus = 0.0, bool i_one_way = false, int *i_flags_ptr = NULL, int i_logger = -1) : bases::boundary (&i_data_plus, i_alpha_plus, &i_data_minus, i_alpha_minus, i_flags_ptr, i_logger) {one_way = i_one_way;}

		boundary (double& i_data_plus, double i_alpha_plus = 0.0, int *i_flags_ptr = NULL, int i_logger = -1) : bases::boundary (&i_data_plus, i_alpha_plus, NULL, 0.0, i_flags_ptr, i_logger) {one_way = false;}
	
		virtual ~boundary () {}
	
		/*!*******************************************************************
		 * \copydoc bases::boundary::execute ()
		 *********************************************************************/
		inline virtual void execute () {
			TRACE (logger, "Executing...");
		
			if (!data_minus) {
				*data_plus *= alpha_plus;
			} else if (!data_plus) {
				*data_minus *= alpha_minus;
			} else {
				*data_minus = alpha_plus * *data_plus + alpha_minus * *data_minus;
				if (!one_way) {
					*data_plus = *data_minus;					
				}
			}
		
			TRACE (logger, "executed.")
		}
	private:
		bool one_way;
	};
} /* one_d */


#endif /* end of include guard: BOUNDARY_HPP_CP5DR4CP */
