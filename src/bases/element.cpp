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
		
		explicit_reset ();
		implicit_reset ();
		
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
	
	void element::update () {
		matrix_solver->update ();
	}
		
	void element::update_timestep (double new_timestep) {
		duration += timestep;
		INFO ("TOTAL TIME: " << duration);
		if (new_timestep != timestep) {
			flags &= ~unchanged_timestep;
			flags &= ~factorized;
			INFO ("Updating timestep: " << timestep);
		} else {
			flags |= unchanged_timestep;
		}
		timestep = new_timestep;
		
		TRACE ("Update complete");
	}
	
	void element::run () {
		double t_timestep;
		for (int i = 0; i < inputParams ["timesteps"].asInt; ++i) {
			INFO ("Timestep " << i);
			DEBUG ("FLAGS " << flags);
			calculate ();
			DEBUG ("FLAGS " << flags);
			output ();
			DEBUG ("FLAGS " << flags);
			execute_boundaries ();
			t_timestep = calculate_timestep ();
			messenger_ptr->min (&t_timestep);
			TRACE ("Updating...");
			for (int k = 0; k < 2; ++k) {
				attempt_update ();
			}
			attempt_update ();
			update ();
			update_timestep (t_timestep);
		}
	}
} /* bases */