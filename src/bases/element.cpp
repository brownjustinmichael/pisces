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
	
/*
	void element::send () {
		for (int i = 0; i < (int) boundary_bools.size (); ++i) {
			send (i);
		}
	}
	
	void element::send (int edge) {
		TRACE (logger, "Sending information...");
		
		for (iterator iter = begin (); iter != end (); ++iter) {
			send (edge, *iter);
		}
		
		TRACE (logger, "Information sent.");
	}
	
	void element::recv () {
		for (int i = 0; i < (int) boundary_bools.size (); ++i) {
			recv (i);
		}
	}

	void element::recv (int edge) {
		TRACE (logger, "Receiving information...");
		
		for (iterator iter = begin (); iter != end (); ++iter) {
			recv (edge, *iter);
		}
		
		TRACE (logger, "Information received.");
	}*/

		
	void element::output () {
		TRACE (logger, "Writing to file...");
	
		// Output in normal space
		if (normal_stream) {
			normal_stream->to_file ();
		}
	}
	
	void element::attempt_update () {
		TRACE (logger, "Updating...");

		if (matrix_solver) {
			matrix_solver->execute ();
		} else {
			WARN (logger, "No matrix solver defined. It is likely the element was not set up correctly.")
		}
		
		TRACE (logger, "Updated.");
	}
	
	void element::calculate_bounds () {
		matrix_solver->calculate_bounds ();
	}
	
	void element::send_bounds () {
		matrix_solver->send_bounds ();
	}
	
	void element::recv_bounds () {
		matrix_solver->recv_bounds ();
	}
	
	void element::calculate_error () {
		matrix_solver->calculate_error ();
	}
	
	void element::send_error () {
		matrix_solver->send_error ();
	}
	
	void element::recv_error () {
		matrix_solver->recv_error ();
	}
	
	void element::update () {
		matrix_solver->update ();
	}
		
	void element::update_timestep (double new_timestep) {
		duration += timestep;
		INFO (logger, "TOTAL TIME: " << duration);
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