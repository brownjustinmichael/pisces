/*!**********************************************************************
 * \file boussinesq.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "boussinesq.hpp"

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

namespace data
{
	template <class datatype>
	thermo_compositional_data <datatype>::thermo_compositional_data (grids::axis *i_axis_n, grids::axis *i_axis_m, int id, int n_elements, io::parameters& i_params) : implemented_data <datatype> (i_axis_n, i_axis_m, i_params, id, i_params.get <std::string> ("dump.file"), i_params.get <std::string> ("root") + i_params.get <std::string> ("dump.directory"), i_params.get <int> ("dump.every")) {
		TRACE ("Initializing...");

		initialize ("pressure", corrector);
		initialize ("composition");
		initialize ("temperature");
		initialize ("x_velocity");
		initialize ("z_velocity");
		initialize ("density");

		// Set up the data from the input file in params
		TRACE ("Setting up from input...");
		i_params ["input.directory"] = i_params ["root"].as <std::string> () + "/" + i_params ["input.directory"].as <std::string> ();
		this->template setup_from <formats::netcdf> (i_params ["input"]);
		
		i_params ["output.directory"] = i_params ["root"].as <std::string> () + "/" + i_params ["output.directory"].as <std::string> ();

		TRACE ("Setting up cartesian output...");
		this->template setup_output_from <formats::netcdf> (i_params ["output.cart"]);

		TRACE ("Setting up transformed output...");
		this->template setup_output_from <formats::netcdf> (i_params ["output.trans"], real_spectral);

		TRACE ("Setting up stat output...");
		std::shared_ptr <io::output> stat_stream =
		this->template setup_output_from <formats::netcdf> (i_params ["output.stat"], real_real, no_variables);

		if (!stat_stream) return;

		for (typename data <datatype>::iterator iter = this->begin (); iter != this->end (); ++iter) {
			// For each data variable, output z_flux, average derivative across the center, average and max
			std::string variable = *iter;
			stat_stream->template append <datatype> ("max_" + variable, this->output_max (variable), formats::scalar);
			stat_stream->template append <datatype> ("avg_" + variable, this->output_avg (variable), formats::scalar);
			stat_stream->template append <datatype> ("deriv_" + variable, this->output_deriv (variable), formats::scalar);
			stat_stream->template append <datatype> ("flux_" + variable, this->output_flux (variable, "z_velocity"), formats::scalar);
		}
	}
	
	template class thermo_compositional_data <double>;
} /* data */

namespace pisces
{
	using namespace plans;
	using namespace plans::solvers;
	using namespace boundaries;
	
	template <class datatype>
	boussinesq_element <datatype>::boussinesq_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data <datatype> &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) : 
	implemented_element <datatype> (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags) {
		TRACE ("Initializing...");
		x_ptr = data ("x");
		z_ptr = data ("z");
		x_vel_ptr = data ("x_velocity", real_real);
		z_vel_ptr = data ("z_velocity", real_real);
		
		cfl = i_params ["time.cfl"].as <datatype> ();
		allow = i_params ["time.allow"].as <datatype> ();

		// Add a background temperature gradient of the form
		// -C*Aout * arctan((z-rt)/dout), z < rt
		// -C*Ain * arctan((z-rt)/din), z > rt
		data.initialize ("korre_Ts", uniform_n);
		if (i_params ["equations.temperature.korre_Ts"].IsDefined ()) {
			datatype C = i_params ["equations.temperature.korre_Ts.C"].as <datatype> ();
			datatype Ain = i_params ["equations.temperature.korre_Ts.Ain"].as <datatype> ();
			datatype din = i_params ["equations.temperature.korre_Ts.din"].as <datatype> ();
			datatype rt = i_params ["equations.temperature.korre_Ts.rt"].as <datatype> ();
			datatype Aout = i_params ["equations.temperature.korre_Ts.Aout"].as <datatype> ();
			datatype dout = i_params ["equations.temperature.korre_Ts.dout"].as <datatype> ();

			for (int j = 0; j < m; ++j)
			{
				if (z_ptr [j] > rt) {
					data ["korre_Ts"] [j] = -C * Aout * std::tanh ((z_ptr [j] - rt) / dout);
				} else {
					data ["korre_Ts"] [j] = -C * Ain * std::tanh ((z_ptr [j] - rt) / din);
				}
			}
		}

		bool uniform_diff = true;
		data.initialize ("temperature_diffusion", uniform_n);
		if (i_params ["equations.temperature.korre_diff"].IsDefined ()) {
			uniform_diff = false;
			datatype Prcz = 1.0 / i_params ["equations.temperature.diffusion"].as <datatype> ();
			datatype Prrz = 1.0 / i_params ["equations.temperature.korre_diff.rz_diffusion"].as <datatype> ();

			DEBUG ("VARS ARE " << Prcz << " " << Prrz);

			datatype A = (Prcz * data ["korre_Ts"] [0] + Prrz) / (data ["korre_Ts"] [0] + 1.);
			assert (A < Prcz);
			for (int j = 0; j < m; ++j)
			{
				data ["temperature_diffusion"] [j] = 1. / (A - data ["korre_Ts"] [j] * (Prcz - A));
				DEBUG ("DIFF IS " << data ["temperature_diffusion"] [j]);
			}
		}

		
		// Set up the temperature equation
		if (!(i_params ["equations.temperature.ignore"].IsDefined () && i_params ["equations.temperature.ignore"].as <bool> ())) {
			*split_solver <datatype> (equations ["temperature"], timestep, 
				neumann (i_params ["equations.temperature.bottom.value"].as <datatype> ()), 
				dirichlet (i_params ["equations.temperature.top.value"].as <datatype> ())) 
			+ advec (data ["x_velocity"], data ["z_velocity"])
			+ src (data ["z_velocity"] * data ["korre_Ts"])
			== 
			params ["equations.temperature.sources.z_velocity"] * src (data ["z_velocity"]);

			if (uniform_diff) {
				*equations ["temperature"] == params ["equations.temperature.diffusion"] * diff ();
			} else {
				*equations ["temperature"] == bg_diff (data ["temperature_diffusion"].ptr ());
			}
			// if (i_params ["equations.temperature.linear"].IsDefined ()) *equations ["temperature"] == plans::diffusion::linear <datatype>::factory (i_params ["equations.temperature.linear"].as <datatype> (), 0.0, data ["temperature_diffusion"].ptr (), 10000);
		}

		// Set up the composition equation
		if (i_params ["equations.composition"].IsDefined () && !(i_params ["equations.composition.ignore"].IsDefined () && i_params ["equations.composition.ignore"].as <bool> ())) {
			*split_solver <datatype> (equations ["composition"], timestep, 
				dirichlet (i_params ["equations.composition.bottom.value"].as <datatype> ()), 
				dirichlet (i_params ["equations.composition.top.value"].as <datatype> ())) 
			+ advec (data ["x_velocity"], data ["z_velocity"]) 
			+ params ["equations.composition.sources.z_velocity"] * src (data ["z_velocity"]) 
			== 
			params ["equations.composition.diffusion"] * diff ();
		}

		// Set up the x_velocity equation, note the missing pressure term, which is handled in div
		if (!(i_params ["equations.x_velocity.ignore"].IsDefined () && i_params ["equations.x_velocity.ignore"].as <bool> ())) {
			*split_solver <datatype> (equations ["x_velocity"], timestep, 
				neumann (0.0), 
				neumann (0.0)) 
			+ advec (data ["x_velocity"], data ["z_velocity"]) 
			== 
			params ["equations.velocity.diffusion"] * diff ();
			if (params.get ("equations.x_velocity.ignore_net", false)) data ["x_velocity"].component_flags |= ignore_net;
		}

		// Set up the z_velocity equation, note the missing pressure term, which is handled in div
		if (!(i_params ["equations.z_velocity.ignore"].IsDefined () && i_params ["equations.z_velocity.ignore"].as <bool> ())) {
			*split_solver <datatype> (equations ["z_velocity"], timestep, 
				dirichlet (0.0), 
				dirichlet (0.0)) 
			+ advec (data ["x_velocity"], data ["z_velocity"]) 
			== 
			params ["equations.z_velocity.sources.temperature"] * src (data ["temperature"])
			+ params ["equations.z_velocity.sources.composition"] * src (data ["composition"]) 
			+ params ["equations.velocity.diffusion"] * diff ();
			if (params.get ("equations.z_velocity.ignore_net", false)) data ["z_velocity"].component_flags |= ignore_net;
		}

		// Set up the velocity constraint
		if (!(i_params ["equations.pressure.ignore"].IsDefined () && i_params ["equations.pressure.ignore"].as <bool> ())) {
			*div <datatype> (equations ["pressure"], equations ["x_velocity"], equations ["z_velocity"])
			==
			0.0;
		}
		
	TRACE ("Initialized.");
	}
	
	template <class datatype>
	datatype boussinesq_element <datatype>::calculate_timestep (int i, int j, formats::virtual_file *virtual_file) {
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
				return std::min (std::abs ((x_ptr [(i + 1) * m + j] - x_ptr [(i - 1) * m + j]) / x_vel_ptr [i * m + j]) * cfl, std::abs ((z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / z_vel_ptr [i * m + j]) * cfl);
			}
		}
		return 1.0 / 0.0;
	}
	
	template class boussinesq_element <double>;
} /* pisces */
