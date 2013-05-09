/*!***********************************************************************
 * \file bases/boundary.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BOUNDARY_HPP_70XV2K58
#define BOUNDARY_HPP_70XV2K58

#include <vector>
#include <functional>
#include <map>
#include "../config.hpp"
#include "plan.hpp"

enum boundary_flags {
	linked_0 = 0x20,
	linked_n = 0x40
};

namespace bases
{
	class element;
	
	class boundary : public plan
	{
	public:
		boundary (int i_edge) : plan () {
			edge = i_edge;
		}
		
		virtual ~boundary () {}
		
		virtual void associate (element* element_ptr);
		
		virtual void send () {}
		
		virtual void recv () {}
	
	protected:
		int edge;
		int index;
		int increment;
	};
} /* bases */

#endif /* end of include guard: BOUNDARY_HPP_70XV2K58 */
