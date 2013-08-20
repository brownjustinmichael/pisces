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

/*!*******************************************************************
 * \brief An enumeration containing the flags for the transform class
 *********************************************************************/
enum transform_flags {
	transformed = 0x10
};

namespace bases
{
	class element;
	
	/*!*******************************************************************
	 * \brief An explicit plan that operates a transform
	 *********************************************************************/
	class transform : public explicit_plan
	{
	public:
		/*!*******************************************************************
		 * \copydoc bases::explicit_plan::explicit_plan ()
		 *********************************************************************/
		transform (element* i_element_ptr, int i_n, int i_name_in, int i_name_out = null) : 
		bases::explicit_plan (i_element_ptr, i_n, i_name_in, i_name_out) {
			TRACE ("Instantiated.");
		}
		
		virtual ~transform () {}
		
		/*!*******************************************************************
		 * \copybrief bases::explicit_plan::execute ()
		 *********************************************************************/
		virtual void execute () {
			explicit_plan::execute ();
		}
	};
} /* bases */

#endif /* end of include guard: TRANSFORM_HPP_27JF6IZD */
