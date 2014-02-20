/*!***********************************************************************
 * \file bases/transform.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2014-02-13.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef TRANSFORM_HPP_S8KRHTX3
#define TRANSFORM_HPP_S8KRHTX3
	
#include "grid.hpp"

/*!**********************************************************************
 * \brief A set of flags to be used with setting up transforms
 ************************************************************************/
enum transform_flags {
	forward_horizontal = 0x01,
	forward_vertical = 0x02,
	inverse_horizontal = 0x04,
	inverse_vertical = 0x08,
	ignore_m = 0x10,
	inverse = 0x20,
};

namespace bases
{
	/*!*******************************************************************
	 * \brief A class designed to solve a matrix equation
	 *********************************************************************/
	template <class datatype>
	class master_transform
	{
	public:
		master_transform (int *i_element_flags, int *i_component_flags) :
		element_flags (i_component_flags),
		component_flags (i_component_flags) {}
		
		/*!*******************************************************************
		 * \brief Solve the matrix equation
		 *********************************************************************/
		virtual void transform (int flags) = 0;
		
		virtual void write () = 0;
		
		virtual void read () = 0;
		
		int *element_flags;
		int *component_flags;
	};
} /* bases */

#endif /* end of include guard: TRANSFORM_HPP_S8KRHTX3 */
