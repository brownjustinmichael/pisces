/*!***********************************************************************
 * \file boundary.hpp
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
#include "../config.hpp"
#include "../plan.hpp"

namespace bases
{
   /*!*******************************************************************
    * \brief An abstract plan class to handle boundary conditions
    * 
    * This class contains the necessary coefficients and pointers to 
    * handle boundary conditions. The implementation will depend on the 
    * number of dimensions.
    *********************************************************************/
	class boundary : public plan
	{
	public:
		/*!*******************************************************************
		 * \param i_alpha_plus A double coefficient for the contribution from the positive boundary
		 * \param i_data_plus A pointer to the double first element of the positive boundary
		 * \param i_alpha_minus A double coefficient for the contribution from the negative boudary
		 * \param i_data_minus A pointer to the double first element of the negative boundary
		 *********************************************************************/
		boundary (double i_alpha_plus, double *i_data_plus, double i_alpha_minus, double *i_data_minus) {
			alpha_plus = i_alpha_plus;
			alpha_minus = i_alpha_minus;
			data_plus = i_data_plus;
			data_minus = i_data_minus;
		}
	
		virtual ~boundary () {}
	
		/*!*******************************************************************
		 * @copydoc plan::execute ()
		 *********************************************************************/
		virtual void execute () = 0;

	protected:
		double alpha_plus; //!< A double coefficient for the contribution from the positive boundary
		double alpha_minus; //!< A double coefficient for the contribution from the negative boudary
		double *data_plus; //!< A pointer to the double first element of the positive boundary
		double *data_minus; //!< A pointer to the double first element of the negative boundary
	};
} /* bases */

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
		 * @copydoc boundary::boundary ()
		 *********************************************************************/
		boundary (double i_alpha_plus, double *i_data_plus, double i_alpha_minus = 0, double *i_data_minus = NULL) : bases::boundary (i_alpha_plus, i_data_plus, i_alpha_minus, i_data_minus) {}
	
		virtual ~boundary () {}
	
		/*!*******************************************************************
		 * @copydoc p::execute ()
		 *********************************************************************/
		inline virtual void execute () {
			TRACE ("Executing...");
		
			if (!data_minus) {
				*data_plus *= alpha_plus;
			} else if (!data_plus) {
				*data_minus *= alpha_minus;
			} else {
				DEBUG ("before = " << *data_plus);
				*data_plus = alpha_plus * *data_plus + alpha_minus * *data_minus;
				DEBUG ("after = " << *data_plus);
				*data_minus = *data_plus;
			}
		
			TRACE ("Executed.")
		}
	
		/*!*******************************************************************
		 * \brief Make a unique pointer to a new boundary_1D instance
		 * 
		 * @copydetails boundary_1D::boundary_1D ()
		 *********************************************************************/
		inline static std::unique_ptr<plan> make_unique (double i_alpha_plus, double *i_data_plus, double i_alpha_minus, double *i_data_minus) {
			return std::unique_ptr<plan> (new boundary (i_alpha_plus, i_data_plus, i_alpha_minus, i_data_minus));
		}
	};
} /* one_d */


#endif /* end of include guard: BOUNDARY_HPP_CP5DR4CP */
