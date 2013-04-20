/*!***********************************************************************
 * \file bases/element.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "element.hpp"
#include "../config.hpp"
#include "solver.hpp"

namespace bases
{
	void element::calculate () {		
		TRACE ("Calculating...");
				
		try {
			// Start in grid space
			
			TRACE ("Updating timestep...");
			
			// Testing
			// Should be replaced by a CFL check
			timestep = 0.001;
			
			TRACE ("Executing explicit grid plans...");
			
			for (int i = 0; i < n_explicit_grid_plans; ++i) {
				explicit_grid_plans [i]->execute ();
			}
			
			TRACE ("Transforming to normal space...");
			
			// Switch to normal space
			if (transform_forward) {
				transform_forward->execute ();
			} else {
				WARN ("Transform not defined. It is likely the element was not set up correctly.")
			}
			
			TRACE ("Writing to file...");
			
			// Output in angle space
			if (angle_stream) {
				angle_stream->to_file ();
			}
			
			TRACE ("Executing explicit space plans...");

			for (int i = 0; i < n_explicit_space_plans; ++i) {
				explicit_space_plans [i]->execute ();
			}
			
			if (timestep != previous_timestep) {
				TRACE ("Executing implicit plans...");
				
				flags &= ~factorized;
				for (int i = 0; i < n_implicit_plans; ++i) {
					implicit_plans [i]->execute ();
				}
			}

			previous_timestep = timestep;
		} catch (...) {
			ERROR ("Exception caught, failsafe dump to _dump.dat");
			failsafe_dump->to_file ();
			exit (EXIT_FAILURE);
		}
		
		TRACE ("Calculation complete.");
	}
		
	void element::execute_boundaries () {
		TRACE ("Executing boundaries...");
		
		for (int i = 0; i < n_boundaries; ++i) {
			boundaries [i]->execute ();
		}
		
		TRACE ("Boundaries executed.");
	}
	
	void element::update () {
		TRACE ("Updating...");
		
		if (matrix_solver) {
			matrix_solver->solve ();
		} else {
			WARN ("No matrix solver defined. It is likely the element was not set up correctly.")
		}
		
		TRACE ("Update complete");
	}
} /* bases */