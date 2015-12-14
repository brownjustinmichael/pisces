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

		// Add a background temperature gradient of the form
		// -C*Aout * arctan((z-rt)/dout), z < rt
		// -C*Ain * arctan((z-rt)/din), z > rt
		data.initialize ("korre_Ts", uniform_n);
		if (i_params ["equations.temperature.korre_Ts"].IsDefined ()) {
			double C = i_params ["equations.temperature.korre_Ts.C"].as <double> ();
			double Ain = i_params ["equations.temperature.korre_Ts.Ain"].as <double> ();
			double din = i_params ["equations.temperature.korre_Ts.din"].as <double> ();
			double rt = i_params ["equations.temperature.korre_Ts.rt"].as <double> ();
			double Aout = i_params ["equations.temperature.korre_Ts.Aout"].as <double> ();
			double dout = i_params ["equations.temperature.korre_Ts.dout"].as <double> ();

			for (int j = 0; j < m; ++j)
			{
				if (z_ptr [j] > rt) {
					data ["korre_Ts"] [j] = -C * Aout * std::tanh ((z_ptr [j] - rt) / dout);
				} else {
					data ["korre_Ts"] [j] = -C * Ain * std::tanh ((z_ptr [j] - rt) / din);
				}
			}
		}

		data.initialize ("temperature_diffusion", uniform_n);
		if (i_params ["equations.temperature.korre_diff"].IsDefined ()) {
			double Prcz = 1.0 / i_params ["equations.temperature.diffusion"].as <double> ();
			double Prrz = 1.0 / i_params ["equations.temperature.korre_diff.rz_diffusion"].as <double> ();

			DEBUG ("VARS ARE " << Prcz << " " << Prrz);

			double A = (Prcz * data ["korre_Ts"] [0] + Prrz) / (data ["korre_Ts"] [0] + 1.);
			assert (A < Prcz);
			for (int j = 0; j < m; ++j)
			{
				data ["temperature_diffusion"] [j] = 1. / (A - data ["korre_Ts"] [j] * (Prcz - A));
				DEBUG ("DIFF IS " << data ["temperature_diffusion"] [j]);
			}
		} else {
			for (int j = 0; j < m; ++j)
			{
				data ["temperature_diffusion"] [j] = params ["equations.temperature.diffusion"].as <double> ();
			}
		}

		// Set up the temperature equation
		if (!(i_params ["equations.temperature.ignore"].IsDefined () && i_params ["equations.temperature.ignore"].as <bool> ())) {
			*split_solver (equations ["temperature"], timestep, 
				neumann (i_params ["equations.temperature.bottom.value"].as <double> ()), 
				dirichlet (i_params ["equations.temperature.top.value"].as <double> ())) 
			+ params ["equations.temperature.advection"] * advec (data ["x_velocity"], data ["z_velocity"])
			+ src (data ["z_velocity"] * data ["korre_Ts"])
			== 
			params ["equations.temperature.sources.z_velocity"] * src (data ["z_velocity"])
			+ bg_diff (data ["temperature_diffusion"].ptr ());

			if (i_params ["equations.temperature.linear"].IsDefined ()) *equations ["temperature"] == std::shared_ptr <plans::plan::factory> (new plans::diffusion::linear::factory (i_params ["equations.temperature.linear"].as <double> (), 0.0, data ["composition"], data ["temperature_diffusion"].ptr (), 10000));
		}
		
	TRACE ("Initialized.");
	}
} /* pisces */
