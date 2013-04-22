/*!***********************************************************************
 * \file bases/transform.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef TRANSFORM_HPP_27JF6IZD
#define TRANSFORM_HPP_27JF6IZD

#include "plan.hpp"

namespace bases
{
	/*!*******************************************************************
	 * \brief An explicit plan that operates a transform
	 *********************************************************************/
	class transform : public explicit_plan
	{
	public:
		/*!*******************************************************************
		 * \copydoc bases::explicit_plan::explicit_plan ()
		 *********************************************************************/
		transform (int i_n, double *i_data_in, double *i_data_out = NULL, int *i_flags_ptr = NULL) : bases::explicit_plan (i_n, i_data_in, i_data_out, i_flags_ptr) {}
		
		virtual ~transform () {}
		
		/*!*******************************************************************
		 * \copybrief bases::explicit_plan::execute ()
		 *********************************************************************/
		virtual void execute () = 0;
	};
} /* bases */

#endif /* end of include guard: TRANSFORM_HPP_27JF6IZD */
