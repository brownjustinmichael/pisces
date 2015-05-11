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
					datatype diff_width = 0.2;
					if (terms ["diff_width"].IsDefined ()) {
						diff_width = terms ["diff_width"].as <datatype> ();
					}
					for (int i = 0; i < m; ++i) {
						if (ptr ("z") [i] > diff_width) {
							diffusion [variable] [i] = terms ["bg_diffusion"].as <datatype> ();
						} else if (ptr ("z") [i] < -diff_width) {
							continue;
						} else {
							diffusion [variable] [i] += (terms ["bg_diffusion"].as <datatype> () - terms ["diffusion"].as <datatype> ()) * (ptr ("z") [i] + diff_width) / (2.0 * diff_width);
						}
					}
				}
			
				// If a diffusion value is specified, construct the diffusion plans
				equations [variable]->add_plan (typename diffusion::background_vertical <datatype>::factory (i_params.get <datatype> ("time.alpha"), &diffusion [variable] [0]), pre_plan);
				// equations [variable]->add_plan (typename diffusion::explicit_background_vertical <datatype>::factory (&diffusion [variable] [0]), mid_plan);
				equations [variable]->add_plan (typename diffusion::background_horizontal <datatype>::factory (i_params.get <datatype> ("time.alpha"), &diffusion [variable] [0]), mid_plan);
			
				if (terms ["variable_diffusion"].IsDefined ()) {
					DEBUG ("Implementing...");
					equations [variable]->add_plan (typename diffusion::linear <datatype>::factory (terms ["variable_diffusion.coefficient"].as <datatype> (), terms ["variable_diffusion.min"].IsDefined () ? terms ["variable_diffusion.min"].as <datatype> () : -(terms ["diffusion"].as <datatype> ()), ptr (terms ["variable_diffusion.source"].as <std::string> ()), &diffusion [variable] [0], 100), post_plan);
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
