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
				
		flags &= ~implicit_started;
		flags &= ~explicit_started;
		
		TRACE (logger, "Writing to file...");
		
		// Output in transform space
		if (transform_stream) {
			transform_stream->to_file ();
		}
		
		TRACE (logger, "Executing plans...");
		
		for (int i = 0; i < (int) plans.size (); ++i) {
			plans [i]->execute ();
		}
		
		TRACE (logger, "Calculation complete.");
	}
	
	void element::send () {
		TRACE (logger, "Sending information...");
		
		for (int i = 0; i < (int) boundaries.size (); ++i) {
			boundaries [i]->send ();
		}
		
		TRACE (logger, "Information sent.");
	}
	
	void element::recv () {
		TRACE (logger, "Receiving information...");
		
		for (int i = 0; i < (int) boundaries.size (); ++i) {
			boundaries [i]->recv ();
		}
		
		TRACE (logger, "Information received.");
	}
		
	void element::execute_boundaries () {
		TRACE (logger, "Executing boundaries...");
		
		for (int i = 0; i < (int) boundaries.size (); ++i) {
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
		
		double new_timestep = calculate_timestep ();
		
		if (matrix_solver) {
			matrix_solver->execute ();
		} else {
			WARN (logger, "No matrix solver defined. It is likely the element was not set up correctly.")
		}
		
		if (new_timestep != timestep) {
			flags &= ~unchanged_timestep;
			flags &= ~factorized;
		} else {
			flags |= unchanged_timestep;
		}
		timestep = new_timestep;
		
		TRACE (logger, "Update complete");
	}
} /* bases */