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

enum transform_flags {
	transformed = 0x10
};

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
		transform (int i_n, int i_name_in, int i_name_out = null) : bases::explicit_plan (i_n, i_name_in, i_name_out) {}
		
		virtual ~transform () {}
		
		/*!*******************************************************************
		 * \copybrief bases::explicit_plan::execute ()
		 *********************************************************************/
		virtual void execute () {}
	};
} /* bases */

#endif /* end of include guard: TRANSFORM_HPP_27JF6IZD */
