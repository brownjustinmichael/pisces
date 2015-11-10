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

	template <class datatype>
	pseudo_element <datatype>::pseudo_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data <datatype> &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) : 
	implemented_element <datatype> (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags) {
		x_ptr = data ("x");
		z_ptr = data ("z");
		x_vel_ptr = data ("x_velocity");
		z_vel_ptr = data ("z_velocity");
		
		cfl = i_params ["time.cfl"].as <datatype> ();

		data.initialize ("bg_pressure", uniform_n);
		data.initialize ("mmw");

		for (int j = 0; j < m; ++j)
		{
			if ((*grids [1]) [j] > 0.0) {
				data ["bg_pressure"].ptr () [j] = i_params ["equations.pressure.zero"].as<datatype>() + (*grids [1]) [j] * (i_params ["equations.pressure.top"].as<datatype>() - i_params ["equations.pressure.zero"].as<datatype>()) * 2.;
			} else {
				data ["bg_pressure"].ptr () [j] = i_params ["equations.pressure.zero"].as<datatype>() + (*grids [1]) [j] * (i_params ["equations.pressure.zero"].as<datatype>() - i_params ["equations.pressure.bottom"].as<datatype>()) * 2.;
			}
		}

		data ["mmw"] == 1.0 / (1.0 - data ["composition"] * i_params ["equations.constants.mass_ratio"].as <datatype> ());
		data ["density"] == data ["bg_pressure"] / data ["temperature"] * data ["mmw"];
		// The mass ratio is (m1 - m2) / m1: near 1 => m1 >> m2, near 0 => m1 ~= m2; m1 > m2 by definition

		data.transformers ["mmw"]->update ();
		data.transformers ["density"]->update ();

		*split_solver <datatype> (equations ["composition"], timestep, 
			dirichlet (i_params ["equations.composition.bottom.value"].as <datatype> ()), 
			dirichlet (i_params ["equations.composition.top.value"].as <datatype> ())) 
		+ advec <datatype> (data ["x_velocity"], data ["z_velocity"]) 
		== 
		params ["equations.composition.diffusion"] * density_diff <datatype> (data ["density"]);

		*split_solver <datatype> (equations ["temperature"], timestep, 
			dirichlet (i_params ["equations.temperature.bottom.value"].as <datatype> ()), 
			dirichlet (i_params ["equations.temperature.top.value"].as <datatype> ())) 
		+ advec <datatype> (data ["x_velocity"], data ["z_velocity"]) 
		== 
		params ["equations.temperature.diffusion"] * density_diff <datatype> (data ["density"] / data ["mmw"])
		// + params ["equations.velocity.diffusion"] * heat <datatype> (data ["mmw"], data ["x_velocity"], data ["z_velocity"])
		- diverge <datatype> (data ["temperature"], data ["x_velocity"], data ["z_velocity"])
		;

		// Set up the x_velocity equation, note the missing pressure term, which is handled in div
		*split_solver <datatype> (equations ["x_velocity"], timestep, neumann (0.0), neumann (0.0)) 
		+ advec <datatype> (data ["x_velocity"], data ["z_velocity"]) 
		== 
		params ["equations.velocity.diffusion"] * density_diff <datatype> (data ["density"])
		+ params ["equations.velocity.diffusion"] * horizontal_stress (data ["density"], data ["z_velocity"])
		;

		// Set up the z_velocity equation, note the missing pressure term, which is handled in div
		*split_solver <datatype> (equations ["z_velocity"], timestep, dirichlet (0.0), dirichlet (0.0)) 
		+ advec <datatype> (data ["x_velocity"], data ["z_velocity"]) 
		== 
		- grad_z <datatype> (data ["bg_pressure"])
		- params ["equations.z_velocity.sources.density"] * src (data ["density"], true)
		+ params ["equations.velocity.diffusion"] * density_diff <datatype> (data ["density"])
		+ params ["equations.velocity.diffusion"] * vertical_stress (data ["density"], data ["x_velocity"])
		;
		data ["x_velocity"].component_flags |= plans::solvers::ignore_net;
		data ["z_velocity"].component_flags |= plans::solvers::ignore_net;

		// *div <datatype> (equations ["pressure"], equations ["x_velocity"], equations ["z_velocity"])
		// Set up the velocity constraint
		*pdiv <datatype> (equations ["pressure"], equations ["x_velocity"], equations ["z_velocity"], data ["density"].ptr (real_spectral), data ["bg_pressure"].ptr (), i_params ["equations.constants.gamma"].as <datatype> ())
		==
		params ["equations.temperature.diffusion"] * density_diff <datatype> (data ["density"] / data ["mmw"], 0.0)
		+ params ["equations.velocity.diffusion"] * heat <datatype> (data ["mmw"], data ["x_velocity"], data ["z_velocity"])
		;
		
	TRACE ("Initialized.");
	}

	template <class datatype>
	datatype pseudo_element <datatype>::calculate_timestep (int i, int j, formats::virtual_file *virtual_file) {
		if (!x_vel_ptr || !z_vel_ptr) {
			return 1.0 / 0.0;
		}
		if (virtual_file) {
			if (j == 0 || j == virtual_file->dims ["z"] [1] - 1) {
				return 1.0 / 0.0;
			}
			if (i == 0 || i == virtual_file->dims ["z"] [0] - 1) {
				return std::abs ((virtual_file->index <datatype> ("z", i, j + 1) - virtual_file->index <datatype> ("z", i, j - 1)) / virtual_file->index <datatype> ("z_velocity", i, j)) * cfl;
			} else {
				return std::min (std::abs ((virtual_file->index <datatype> ("x", i + 1, j) - virtual_file->index <datatype> ("x", i - 1, j)) / virtual_file->index <datatype> ("x_velocity", i, j)), std::abs ((virtual_file->index <datatype> ("z", i, j + 1) - virtual_file->index <datatype> ("z", i, j - 1)) / virtual_file->index <datatype> ("z_velocity", i, j))) * cfl;
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

	template class pseudo_element <double>;
} /* pisces */
