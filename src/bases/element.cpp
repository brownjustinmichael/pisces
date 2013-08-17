/*!***********************************************************************
 * \file bases/element.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cassert>
#include "element.hpp"
#include "../config.hpp"
#include "solver.hpp"
#include "transform.hpp"
#include "../utils/utils.hpp"

namespace bases
{
	void element::send_positions () {
		matrix_solver->send_positions ();
	}
	
	void element::run () {
		double t_timestep;
		for (int i = 0; i < inputParams ["timesteps"].asInt; ++i) {
			INFO ("Timestep " << i);
			
			TRACE ("Calculating...");
		
			explicit_reset ();
			implicit_reset ();
		
			// Output in transform space
			if (transform_stream) {
				TRACE ("Writing to file...");
				transform_stream->to_file ();
			}
		
			TRACE ("Executing plans...");
		
			for (int i = 0; i < (int) plans.size (); ++i) {
				plans [i]->execute ();
			}
		
			TRACE ("Calculation complete.");
			
			if (normal_stream) {
				TRACE ("Writing to file...");
				normal_stream->to_file ();
			}

			execute_boundaries ();

			t_timestep = calculate_timestep ();
			messenger_ptr->min (&t_timestep);
			
			TRACE ("Attempting update...");
			for (int k = 0; k < 3; ++k) {
				if (matrix_solver) {
					matrix_solver->execute ();
				} else {
					WARN ("No matrix solver defined. It is likely the element was not set up correctly.");
				}
			}
			TRACE ("Updating...")
			if (matrix_solver) {
				matrix_solver->update ();
			} else {
				WARN ("No matrix solver defined. It is likely the element was not set up correctly.");
			}
			
			duration += timestep;
			INFO ("TOTAL TIME: " << duration);
			if (t_timestep != timestep) {
				flags &= ~unchanged_timestep;
				flags &= ~factorized;
				INFO ("Updating timestep: " << t_timestep);
			} else {
				flags |= unchanged_timestep;
			}
			timestep = t_timestep;
		
			TRACE ("Update complete");
		}
	}
} /* bases */