/*!**********************************************************************
 * \file pseudo_element.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "pseudo_element.hpp"

#include "plans/plan.hpp"
#include "plans/advection.hpp"
#include "plans/diffusion.hpp"
#include "plans/source.hpp"
#include "plans-solvers/solvers.hpp"

namespace pisces
{
	using namespace plans;

	implemented_element::registrar <pseudo_element> pseudo_registrar ("pseudo");

	pseudo_element::pseudo_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) : 
	implemented_element (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags) {
		x_ptr = data ("x");
		z_ptr = data ("z");
		x_vel_ptr = data ("x_velocity");
		z_vel_ptr = data ("z_velocity");
		
		cfl = i_params ["time.cfl"].as <double> ();

		data.initialize ("bg_pressure", uniform_n);
		data.initialize ("mmw");

		for (int j = 0; j < m; ++j)
		{
			if ((*grids [1]) [j] > 0.0) {
				data ["bg_pressure"].ptr () [j] = i_params ["equations.pressure.zero"].as<double>() + (*grids [1]) [j] * (i_params ["equations.pressure.top"].as<double>() - i_params ["equations.pressure.zero"].as<double>()) * 2.;
			} else {
				data ["bg_pressure"].ptr () [j] = i_params ["equations.pressure.zero"].as<double>() + (*grids [1]) [j] * (i_params ["equations.pressure.zero"].as<double>() - i_params ["equations.pressure.bottom"].as<double>()) * 2.;
			}
		}

		data ["mmw"] == 1.0 / (1.0 - data ["composition"] * i_params ["equations.constants.mass_ratio"].as <double> ());
		data ["density"] == data ["bg_pressure"] / data ["temperature"] * data ["mmw"];
		// The mass ratio is (m1 - m2) / m1: near 1 => m1 >> m2, near 0 => m1 ~= m2; m1 > m2 by definition

		data.transformers ["mmw"]->update ();
		data.transformers ["density"]->update ();

		*split_solver (equations ["composition"], timestep, 
			dirichlet (i_params ["equations.composition.bottom.value"].as <double> ()), 
			dirichlet (i_params ["equations.composition.top.value"].as <double> ())) 
		+ advec (data ["x_velocity"], data ["z_velocity"]) 
		== 
		params ["equations.composition.diffusion"] * density_diff (data ["density"]);

		*split_solver (equations ["temperature"], timestep, 
			dirichlet (i_params ["equations.temperature.bottom.value"].as <double> ()), 
			dirichlet (i_params ["equations.temperature.top.value"].as <double> ())) 
		+ advec (data ["x_velocity"], data ["z_velocity"]) 
		== 
		params ["equations.temperature.diffusion"] * density_diff (data ["density"] / data ["mmw"])
		// + params ["equations.velocity.diffusion"] * heat <double> (data ["mmw"], data ["x_velocity"], data ["z_velocity"])
		- diverge (data ["temperature"], data ["x_velocity"], data ["z_velocity"])
		;

		// Set up the x_velocity equation, note the missing pressure term, which is handled in div
		*split_solver (equations ["x_velocity"], timestep, neumann (0.0), neumann (0.0)) 
		+ advec (data ["x_velocity"], data ["z_velocity"]) 
		== 
		params ["equations.velocity.diffusion"] * density_diff (data ["density"])
		+ params ["equations.velocity.diffusion"] * horizontal_stress (data ["density"], data ["z_velocity"])
		;

		// Set up the z_velocity equation, note the missing pressure term, which is handled in div
		*split_solver (equations ["z_velocity"], timestep, dirichlet (0.0), dirichlet (0.0)) 
		+ advec (data ["x_velocity"], data ["z_velocity"]) 
		== 
		- grad_z (data ["bg_pressure"])
		- params ["equations.z_velocity.sources.density"] * src (data ["density"], true)
		+ params ["equations.velocity.diffusion"] * density_diff (data ["density"])
		+ params ["equations.velocity.diffusion"] * vertical_stress (data ["density"], data ["x_velocity"])
		;
		data ["x_velocity"].component_flags |= plans::solvers::ignore_net;
		data ["z_velocity"].component_flags |= plans::solvers::ignore_net;

		// *div <double> (equations ["pressure"], equations ["x_velocity"], equations ["z_velocity"])
		// Set up the velocity constraint
		*pdiv <double> (equations ["pressure"], equations ["x_velocity"], equations ["z_velocity"], data ["density"].ptr (real_spectral), data ["bg_pressure"].ptr (), i_params ["equations.constants.gamma"].as <double> ())
		==
		params ["equations.temperature.diffusion"] * density_diff (data ["density"] / data ["mmw"], 0.0)
		+ params ["equations.velocity.diffusion"] * heat (data ["mmw"], data ["x_velocity"], data ["z_velocity"])
		;
		
	TRACE ("Initialized.");
	}

	double pseudo_element::calculate_timestep (int i, int j, formats::virtual_file *virtual_file) {
		if (!x_vel_ptr || !z_vel_ptr) {
			return 1.0 / 0.0;
		}
		if (virtual_file) {
			if (j == 0 || j == virtual_file->dims ["z"] [1] - 1) {
				return 1.0 / 0.0;
			}
			if (i == 0 || i == virtual_file->dims ["z"] [0] - 1) {
				return std::abs ((virtual_file->index <double> ("z", i, j + 1) - virtual_file->index <double> ("z", i, j - 1)) / virtual_file->index <double> ("z_velocity", i, j)) * cfl;
			} else {
				return std::min (std::abs ((virtual_file->index <double> ("x", i + 1, j) - virtual_file->index <double> ("x", i - 1, j)) / virtual_file->index <double> ("x_velocity", i, j)), std::abs ((virtual_file->index <double> ("z", i, j + 1) - virtual_file->index <double> ("z", i, j - 1)) / virtual_file->index <double> ("z_velocity", i, j))) * cfl;
			}
		} else {
			if (j == 0 || j == m - 1) {
				return 1.0 / 0.0;
			}
			if ((i == 0 || i == n - 1)) {
				return std::abs ((z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / z_vel_ptr [i * m + j]) * cfl;
			} else {
				return std::min (std::abs ((x_ptr [(i + 1) * m + j] - x_ptr [(i - 1) * m + j]) / x_vel_ptr [i * m + j]), std::abs ((z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / z_vel_ptr [i * m + j])) * cfl;
			}
		}
		return 1.0 / 0.0;
	}
} /* pisces */
