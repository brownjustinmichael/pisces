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

namespace bases
{
	void element::calculate () {		
		TRACE (logger, "Calculating...");
				
		try {
			// Start in grid space
			
			TRACE (logger, "Writing to file...");
			
			// Output in transform space
			if (transform_stream) {
				transform_stream->to_file ();
			}
			
			TRACE (logger, "Updating timestep...");
			
			// Testing
			// Should be replaced by a CFL check
			timestep = 0.001;
			
			TRACE (logger, "Executing explicit grid plans...");
			
			for (int i = 0; i < n_explicit_grid_plans; ++i) {
				explicit_grid_plans [i]->execute ();
			}
			
			TRACE (logger, "Transforming to normal space...");
			
			// Switch to normal space
			if (transform_forward) {
				transform_forward->execute ();
			} else {
				WARN (logger, "Transform not defined. It is likely the element was not set up correctly.")
			}
			
			TRACE (logger, "Executing explicit space plans...");

			for (int i = 0; i < n_explicit_space_plans; ++i) {
				explicit_space_plans [i]->execute ();
			}
			
			if (timestep != previous_timestep) {
				TRACE (logger, "Executing implicit plans...");
				
				flags &= ~factorized;
				for (int i = 0; i < n_implicit_plans; ++i) {
					implicit_plans [i]->execute ();
				}
			}
			
			previous_timestep = timestep;
		} catch (...) {
			ERROR (logger, "Exception caught, failsafe dump to _dump.dat");
			failsafe_dump->to_file ();
			exit (EXIT_FAILURE);
		}
		
		TRACE (logger, "Calculation complete.");
	}
		
	void element::execute_boundaries () {
		TRACE (logger, "Executing boundaries...");
		
		for (int i = 0; i < n_boundaries; ++i) {
			boundaries [i]->execute ();
		}
		
		TRACE (logger, "Writing to file...");
		
		// Output in normal space
		if (normal_stream) {
			normal_stream->to_file ();
		}
		
		TRACE (logger, "Boundaries executed.");
	}
	
	void element::update () {
		TRACE (logger, "Updating...");
		
		// assert (!(flags & transformed));
		
		if (matrix_solver) {
			matrix_solver->solve ();
		} else {
			WARN (logger, "No matrix solver defined. It is likely the element was not set up correctly.")
		}
		
		TRACE (logger, "Update complete");
	}
} /* bases */