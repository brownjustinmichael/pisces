/*!***********************************************************************
 * \file bases/element.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "element.hpp"

#include <cassert>
#include "../config.hpp"
#include "solver.hpp"
#include "../utils/utils.hpp"

namespace bases
{
	template <class datatype>
	void element <datatype>::run () {
		TRACE ("Running...");
		
		for (int j = 0; j < params.timesteps; ++j) {
			INFO ("Timestep " << j << " of " << params.timesteps);
			
			TRACE ("Calculating...");
		
			explicit_reset ();
		
			TRACE ("Executing plans...");
		
			for (int i = 0; i < (int) pre_transform_plans.size (); ++i) {
				pre_transform_plans [i]->execute (flags);
			}
			
			transform_vertical_inverse ();
			
			for (int i = 0; i < (int) mid_transform_plans.size (); ++i) {
				mid_transform_plans [i]->execute (flags);
			}
			
			transform_horizontal_inverse ();

			factorize ();
			
			for (int i = 0; i < (int) post_transform_plans.size (); ++i) {
				post_transform_plans [i]->execute (flags);
			}
		
			TRACE ("Calculation complete.");

			
			if (normal_stream) {
				TRACE ("Writing to file...");
				normal_stream->to_file ();
			}
			
			TRACE ("Updating...");
			
			transform_horizontal_forward ();
			
			solve ();
			
			// Output in transform space
			if (transform_stream) {
				TRACE ("Writing to file...");
				transform_stream->to_file ();
			}
		
			TRACE ("Update complete");
		}
	}
	
	template class element <double>;
} /* bases */