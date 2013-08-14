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
	void element::calculate () {		
		TRACE ("Calculating...");
				
		flags &= ~implicit_started;
		flags &= ~explicit_started;
		
		TRACE ("Writing to file...");
		
		// Output in transform space
		if (transform_stream) {
			transform_stream->to_file ();
		}
		
		TRACE ("Executing plans...");
		
		for (int i = 0; i < (int) plans.size (); ++i) {
			plans [i]->execute ();
		}
		
		TRACE ("Calculation complete.");
	}
	
	void element::output () {
		TRACE ("Writing to file...");
	
		// Output in normal space
		if (normal_stream) {
			normal_stream->to_file ();
		}
	}
	
	void element::attempt_update () {
		TRACE ("Updating...");

		if (matrix_solver) {
			matrix_solver->execute ();
		} else {
			WARN ("No matrix solver defined. It is likely the element was not set up correctly.");
		}
		
		TRACE ("Updated.");
	}
	
	void element::send_positions () {
		matrix_solver->send_positions ();
	}
	
	void element::recv_positions () {
		matrix_solver->recv_positions ();
	}
	
	void element::update () {
		matrix_solver->update ();
	}
		
	void element::update_timestep (double new_timestep) {
		duration += timestep;
		INFO ("TOTAL TIME: " << duration);
		if (new_timestep != timestep) {
			flags &= ~unchanged_timestep;
			flags &= ~factorized;
		} else {
			flags |= unchanged_timestep;
		}
		timestep = new_timestep;
		
		TRACE ("Update complete");
	}
} /* bases */