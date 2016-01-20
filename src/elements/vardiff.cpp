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
#include "plans/advection.hpp"
#include "plans-solvers/boundaries/implemented_boundary.hpp"
#include "plans-solvers/solvers.hpp"


namespace pisces
{
	using namespace plans;
	using namespace plans::solvers;

	implemented_element::registrar <vardiff_element> vardiff_registrar ("vardiff");
	
	vardiff_element::vardiff_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) : 
	boussinesq_element (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags) {
		TRACE ("Initializing...");
		
		data.initialize ("temperature_diffusion", uniform_n);
		for (int j = 0; j < m; ++j)
		{
			data ["temperature_diffusion"] [j] = params ["equations.temperature.diffusion"].as <double> ();
		}

		double phi = i_params ["equations.temperature.tanh.coeff"].as <double> ();
		double stiffness = i_params ["equations.temperature.stiffness"].as <double> ();
		double chi = 2. * phi / (1.0 + phi) / (stiffness * (1.0 + phi) / (1.0 - phi) + 1.0);
		WARN("chi = " << chi);

		if (!(i_params ["equations.temperature.ignore"].IsDefined () && i_params ["equations.temperature.ignore"].as <bool> ())) {
			*split_solver (equations ["temperature"], timestep, 
				dirichlet (i_params ["equations.temperature.bottom.value"].as <double> ()), 
				neumann (-1.0 / chi)) 
			+ params ["equations.temperature.advection"] * advec (data ["x_velocity"], data ["z_velocity"])
			+ (1.0 / chi - 1.0) * src (data ["z_velocity"])
			== 
			params ["equations.temperature.sources.z_velocity"] * src (data ["z_velocity"])
			+ bg_diff (data ["temperature_diffusion"].ptr ());

			if (i_params ["equations.temperature.linear"].IsDefined ()) *equations ["temperature"] == std::shared_ptr <plans::plan::factory> (new plans::diffusion::linear::factory (i_params ["equations.temperature.linear"].as <double> () * i_params ["equations.temperature.diffusion"].as <double> (), -1.0, data ["composition"], data ["temperature_diffusion"].ptr (), 10));
			if (i_params ["equations.temperature.tanh"].IsDefined ()) *equations ["temperature"] == std::shared_ptr <plans::plan::factory> (new plans::diffusion::tanh::factory (i_params ["equations.temperature.diffusion"].as <double> () * i_params ["equations.temperature.tanh.coeff"].as <double> (), i_params ["equations.temperature.tanh.length"].as <double> (), i_params ["equations.temperature.tanh.zero"].as <double> (), -1.0, data ["composition"], data ["temperature_diffusion"].ptr (), 10));

		}
		TRACE ("Initialized.");
	}
} /* pisces */
