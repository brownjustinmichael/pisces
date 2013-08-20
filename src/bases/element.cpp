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
	void element::run () {
		double t_timestep;
		implicit_reset ();

		for (std::shared_ptr <plan> i_plan : implicit_plans) {
			i_plan->execute ();
		}
		
		for (int j = 0; j < inputParams ["timesteps"].asInt; ++j) {
			INFO ("Timestep " << j);
			
			TRACE ("Calculating...");
		
			explicit_reset ();
		
			// Output in transform space
			if (transform_stream) {
				TRACE ("Writing to file...");
				transform_stream->to_file ();
			}
		
			TRACE ("Executing plans...");
		
			for (std::shared_ptr <plan> i_plan : plans) {
				i_plan->execute ();
			}
		
			TRACE ("Calculation complete.");
			
			if (normal_stream) {
				TRACE ("Writing to file...");
				normal_stream->to_file ();
			}

			execute_boundaries ();

			t_timestep = calculate_timestep ();
			messenger_ptr->min (&t_timestep);
			
			TRACE ("Updating...")
			for (std::shared_ptr <solver> i_solver : solvers) {
				i_solver->execute ();
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