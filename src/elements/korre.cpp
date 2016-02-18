/*!**********************************************************************
 * \file korre.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "korre.hpp"

#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <memory>

#include "io/formats/netcdf.hpp"
#include "io/formats/ascii.hpp"

#include "plans/plan.hpp"
#include "plans/advection.hpp"
#include "plans/diffusion.hpp"
#include "plans/source.hpp"
#include "plans-solvers/boundaries/implemented_boundary.hpp"
#include "plans-solvers/solvers.hpp"
#include "io/functors/div.hpp"

namespace pisces
{
	using namespace plans;
	using namespace plans::solvers;
	using namespace boundaries;

	implemented_element::registrar <korre_element> korre_registrar ("korre");
	
	korre_element::korre_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) : 
	boussinesq_element (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags) {
		TRACE ("Initializing...");

		data.initialize ("korre_Ts", uniform_n);
		// data.initialize ("temperature_diffusion", uniform_n);

		double stiffness = i_params ["equations.temperature.stiffness"].as <double> ();
		double diff = i_params ["equations.temperature.diffusion"].as <double> ();
		double chi = 1. / (1. - i_params.get <double> ("equations.temperature.sources.z_velocity"));

		if (i_params ["equations.temperature.tanh"].IsDefined ()) {

			double stiffness = i_params ["equations.temperature.stiffness"].as <double> ();
			double zero = i_params ["equations.temperature.tanh.zero"].as <double> ();
			double length = i_params ["equations.temperature.tanh.length"].as <double> ();

			for (int j = 0; j < m; ++j)
			{
				data ["korre_Ts"] [j] = (1.0 + stiffness) / 2.0 * std::tanh ((z_ptr [j] - zero) / length) + (1.0 - stiffness) / 2.0;
			}
			for (int j = 0; j < m; ++j)
			{
				data ["temperature_diffusion"] [j] = 2.0 * diff / (1.0 - chi + sqrt((chi - 1.0) * (chi - 1.0) + 4.0 * chi * data ["korre_Ts"] [j]));
				DEBUG ("DIFF IS " << data ["temperature_diffusion"] [j]);
			}
		} else {
			for (int j = 0; j < m; ++j)
			{
				data ["temperature_diffusion"] [j] = diff;
			}
		}

		// Set up the temperature equation
		if (!(i_params ["equations.temperature.ignore"].IsDefined () && i_params ["equations.temperature.ignore"].as <bool> ())) {
			*split_solver (equations ["temperature"], timestep, 
				dirichlet (i_params ["equations.temperature.bottom.value"].as <double> ()), 
				neumann (i_params ["equations.temperature.top.value"].as <double> ())) 
			+ params ["equations.temperature.advection"] * advec (data ["x_velocity"], data ["z_velocity"])
			// + src (data ["z_velocity"] * data ["korre_Ts"])
			== 
			params ["equations.temperature.sources.z_velocity"] * src (data ["z_velocity"])
			+ bg_diff (data ["temperature_diffusion"].ptr ());
		}
		
	TRACE ("Initialized.");
	}
} /* pisces */
