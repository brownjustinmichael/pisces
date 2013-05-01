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
#include "element.hpp"

enum boundary_flags {
	fixed_0 = 0x20,
	fixed_n = 0x21
}

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
		boundary (int *i_flags_ptr = NULL, int i_logger = -1) : plan (i_flags_ptr, i_logger) {
		}
	
		virtual ~boundary () {}
	
		/*!*******************************************************************
		 * \copydoc plan::execute ()
		 *********************************************************************/
		virtual void execute () {
			
		}

	protected:

	};
} /* bases */

#endif /* end of include guard: BOUNDARY_HPP_70XV2K58 */
