/*!***********************************************************************
 * \file bases/boundary.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BOUNDARY_HPP_70XV2K58
#define BOUNDARY_HPP_70XV2K58

#include "plan.hpp"

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
		boundary (double i_alpha_plus, double *i_data_plus, double i_alpha_minus, double *i_data_minus, int *i_flags_ptr = NULL, int i_logger = -1) : plan (i_flags_ptr, i_logger) {
			alpha_plus = i_alpha_plus;
			alpha_minus = i_alpha_minus;
			data_plus = i_data_plus;
			data_minus = i_data_minus;
		}
	
		virtual ~boundary () {}
	
		/*!*******************************************************************
		 * \copydoc plan::execute ()
		 *********************************************************************/
		virtual void execute () = 0;

	protected:
		double alpha_plus; //!< A double coefficient for the contribution from the positive boundary
		double alpha_minus; //!< A double coefficient for the contribution from the negative boudary
		double *data_plus; //!< A pointer to the double first element of the positive boundary
		double *data_minus; //!< A pointer to the double first element of the negative boundary
	};
} /* bases */

#endif /* end of include guard: BOUNDARY_HPP_70XV2K58 */
