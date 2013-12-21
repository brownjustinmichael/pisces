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
		
		for (int j = 0; j <= params.timesteps; ++j) {
			INFO ("Timestep " << j << " of " << params.timesteps);
			
			TRACE ("Calculating...");
		
			explicit_reset ();
		
			TRACE ("Executing plans...");
		
			for (iterator iter = begin (); iter != end (); iter++) {
				iter->second->execute_pre_plans (flags [iter->first]);
			}
			
			transform (inverse_vertical);
			
			for (iterator iter = begin (); iter != end (); iter++) {
				iter->second->execute_mid_plans (flags [iter->first]);
			}
			
			transform (inverse_horizontal);

			factorize ();
			
			for (iterator iter = begin (); iter != end (); iter++) {
				iter->second->execute_post_plans (flags [iter->first]);
			}
		
			TRACE ("Calculation complete.");

			
			if (normal_stream) {
				TRACE ("Writing to file...");
				normal_stream->to_file ();
			}
			
			TRACE ("Updating...");

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