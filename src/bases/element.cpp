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
		
		for (int i = 0; i < (int) plans.size (); ++i) {
			plans [i]->execute ();
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
		
		if (matrix_solver) {
			matrix_solver->execute ();
		} else {
			WARN (logger, "No matrix solver defined. It is likely the element was not set up correctly.")
		}
	
		TRACE (logger, "Update complete");
	}
} /* bases */