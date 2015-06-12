/*!**********************************************************************
 * \file vardiff.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "vardiff.hpp"

#include "plans/plan.hpp"
#include "plans/diffusion.hpp"

namespace pisces
{
	using namespace plans;
	using namespace plans::solvers;
	
	template <class datatype>
	vardiff_element <datatype>::vardiff_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data <datatype> &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) : 
	boussinesq_element <datatype> (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags, false) {
		TRACE ("Initializing...");
		
		std::string variable;
		for (YAML::const_iterator iter = i_params ["equations"].begin (); iter != i_params ["equations"].end (); ++iter) {
			variable = iter->first.as <std::string> ();
			INFO ("Setting up " << variable);
			io::parameters::alias terms (i_params, "equations." + variable);
			if (terms ["ignore"].IsDefined () && terms ["ignore"].as <bool> ()) {
				continue;
			}
			
			int m = i_axis_m.get_n ();
			diffusion [variable].resize (m);
			
			if (terms ["diffusion"].IsDefined ()) {
				for (int i = 0; i < m; ++i) {
					diffusion [variable] [i] = terms ["diffusion"].as <datatype> ();
				}
				
				if (terms ["bg_diffusion"].IsDefined ()) {
					datatype lo_diffusion = terms ["diffusion"].as <datatype> ();
					datatype hi_diffusion = terms ["bg_diffusion"].as <datatype> ();
					datatype diff_width = 0.2;
					if (terms ["diff_width"].IsDefined ()) {
						diff_width = terms ["diff_width"].as <datatype> ();
					}
					for (int i = 0; i < m; ++i) {
						diffusion [variable] [i] = exp (atan (ptr ("z") [i] / diff_width) * (log (hi_diffusion) - log(lo_diffusion)) / 3.14159 + log (lo_diffusion) + (log (hi_diffusion) - log (lo_diffusion)) / 2.0);
						DEBUG ("Diffusion: " << diffusion [variable] [i]);
					}
				}
			
				// If a diffusion value is specified, construct the diffusion plans
				equations [variable]->add_plan (typename diffusion::background_vertical <datatype>::factory (i_params.get <datatype> ("time.alpha"), &diffusion [variable] [0]));
				equations [variable]->add_plan (typename diffusion::background_horizontal <datatype>::factory (i_params.get <datatype> ("time.alpha"), &diffusion [variable] [0]));
			
				if (terms ["variable_diffusion"].IsDefined ()) {
					DEBUG ("Implementing...");
					equations [variable]->add_plan (typename diffusion::linear <datatype>::factory (terms ["variable_diffusion.coefficient"].as <datatype> (), terms ["variable_diffusion.min"].IsDefined () ? terms ["variable_diffusion.min"].as <datatype> () : -(terms ["diffusion"].as <datatype> ()), ptr (terms ["variable_diffusion.source"].as <std::string> ()), &diffusion [variable] [0], 100));
				}
			} else if (terms ["variable_diffusion"].IsDefined ()) {
				FATAL ("Can't use variable diffusion without base diffusion yet.");
				throw 0;
			}

		}
		TRACE ("Initialized.");
	}
	
	template class vardiff_element <double>;
} /* pisces */
