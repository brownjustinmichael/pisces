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
#include "io/functors/div.hpp"

namespace data
{
	thermo_compositional_data::thermo_compositional_data (grids::axis *i_axis_n, grids::axis *i_axis_m, int id, int n_elements, io::parameters& i_params) : implemented_data (i_axis_n, i_axis_m, i_params, id, i_params.get <std::string> ("dump.file"), i_params.get <std::string> ("root") + i_params.get <std::string> ("dump.directory"), i_params.get <int> ("dump.every")) {
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
		std::shared_ptr <io::output> stream = this->template setup_output_from <formats::netcdf> (i_params ["output.cart"]);
		// stream->template append <double> ("div_u", std::shared_ptr <typename functors::div_functor <double>> (new functors::div_functor <double> ((*this) ["x"].ptr (), (*this) ["z"].ptr (), (*this) ["x_velocity"].ptr (), (*this) ["z_velocity"].ptr (), n, m)));

		TRACE ("Setting up transformed output...");
		stream = this->template setup_output_from <formats::netcdf> (i_params ["output.trans"], real_spectral);
		// stream->template append <double> ("div_u", std::shared_ptr <typename functors::transform_div_functor <double>> (new functors::transform_div_functor <double> ((*this) ["x"].ptr (), (*this) ["z"].ptr (), (*this) ["x_velocity"].ptr (), (*this) ["z_velocity"].ptr (), n, m)));

		TRACE ("Setting up stat output...");
		std::shared_ptr <io::output> stat_stream =
		this->template setup_output_from <formats::netcdf> (i_params ["output.stat"], real_real, no_variables);

		if (!stat_stream) return;

		for (typename data::iterator iter = this->begin (); iter != this->end (); ++iter) {
			// For each data variable, output z_flux, average derivative across the center, average and max
			std::string variable = *iter;
			stat_stream->template append <double> ("max_" + variable, this->output_max (variable), formats::scalar);
			stat_stream->template append <double> ("avg_" + variable, this->output_avg (variable), formats::scalar);
			stat_stream->template append <double> ("deriv_" + variable, this->output_deriv (variable), formats::scalar);
			stat_stream->template append <double> ("flux_" + variable, this->output_flux (variable, "z_velocity"), formats::scalar);
		}
	}
} /* data */

namespace pisces
{
	using namespace plans;
	using namespace plans::solvers;
	using namespace boundaries;
	
	boussinesq_element::boussinesq_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) : 
	implemented_element (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags) {
		TRACE ("Initializing...");
		x_ptr = data ("x");
		z_ptr = data ("z");
		x_vel_ptr = data ("x_velocity", real_real);
		z_vel_ptr = data ("z_velocity", real_real);
		
		cfl = i_params ["time.cfl"].as <double> ();
		allow = i_params ["time.allow"].as <double> ();

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

		bool uniform_diff = true;
		// data.initialize ("temperature_diffusion", uniform_n);
		// if (i_params ["equations.temperature.korre_diff"].IsDefined ()) {
		// 	uniform_diff = false;
		// 	double Prcz = 1.0 / i_params ["equations.temperature.diffusion"].as <double> ();
		// 	double Prrz = 1.0 / i_params ["equations.temperature.korre_diff.rz_diffusion"].as <double> ();

		// 	DEBUG ("VARS ARE " << Prcz << " " << Prrz);

		// 	double A = (Prcz * data ["korre_Ts"] [0] + Prrz) / (data ["korre_Ts"] [0] + 1.);
		// 	assert (A < Prcz);
		// 	for (int j = 0; j < m; ++j)
		// 	{
		// 		data ["temperature_diffusion"] [j] = 1. / (A - data ["korre_Ts"] [j] * (Prcz - A));
		// 		DEBUG ("DIFF IS " << data ["temperature_diffusion"] [j]);
		// 	}
		// }

		
		// Set up the temperature equation
		if (!(i_params ["equations.temperature.ignore"].IsDefined () && i_params ["equations.temperature.ignore"].as <bool> ())) {
			*split_solver (equations ["temperature"], timestep, 
				dirichlet (i_params ["equations.temperature.bottom.value"].as <double> ()), 
				dirichlet (i_params ["equations.temperature.top.value"].as <double> ())) 
			+ params ["equations.temperature.advection"] * advec (data ["x_velocity"], data ["z_velocity"])
			// + src (data ["z_velocity"] * data ["korre_Ts"])
			== 
			params ["equations.temperature.sources.z_velocity"] * src (data ["z_velocity"]);

			if (uniform_diff) {
				*equations ["temperature"] == params ["equations.temperature.diffusion"] * diff ();
			} else {
				*equations ["temperature"] == bg_diff (data ["temperature_diffusion"].ptr ());
			}
			// if (i_params ["equations.temperature.linear"].IsDefined ()) *equations ["temperature"] == plans::diffusion::linear <double>::factory (i_params ["equations.temperature.linear"].as <double> (), 0.0, data ["temperature_diffusion"].ptr (), 10000);
		}

		// Set up the composition equation
		if (i_params ["equations.composition"].IsDefined () && !(i_params ["equations.composition.ignore"].IsDefined () && i_params ["equations.composition.ignore"].as <bool> ())) {
			*split_solver (equations ["composition"], timestep, 
				dirichlet (i_params ["equations.composition.bottom.value"].as <double> ()), 
				dirichlet (i_params ["equations.composition.top.value"].as <double> ())) 
			+ params ["equations.composition.advection"] * advec (data ["x_velocity"], data ["z_velocity"]) 
			+ params ["equations.composition.sources.z_velocity"] * src (data ["z_velocity"]) 
			== 
			params ["equations.composition.diffusion"] * diff ();
		}

		// Set up the x_velocity equation, note the missing pressure term, which is handled in div
		if (!(i_params ["equations.x_velocity.ignore"].IsDefined () && i_params ["equations.x_velocity.ignore"].as <bool> ())) {
			*split_solver (equations ["x_velocity"], timestep, 
				neumann (0.0), 
				neumann (0.0)) 
			+ params ["equations.velocity.advection"] * advec (data ["x_velocity"], data ["z_velocity"]) 
			== 
			params ["equations.velocity.diffusion"] * diff ();
			if (params.get ("equations.x_velocity.ignore_net", false)) data ["x_velocity"].component_flags |= ignore_net;
		}

		// Set up the z_velocity equation, note the missing pressure term, which is handled in div
		if (!(i_params ["equations.z_velocity.ignore"].IsDefined () && i_params ["equations.z_velocity.ignore"].as <bool> ())) {
			*split_solver (equations ["z_velocity"], timestep, 
				dirichlet (0.0), 
				dirichlet (0.0)) 
			+ params ["equations.velocity.advection"] * advec (data ["x_velocity"], data ["z_velocity"]) 
			== 
			params ["equations.z_velocity.sources.temperature"] * src (data ["temperature"])
			+ params ["equations.z_velocity.sources.composition"] * src (data ["composition"]) 
			+ params ["equations.velocity.diffusion"] * diff ();
			if (params.get ("equations.z_velocity.ignore_net", false)) data ["z_velocity"].component_flags |= ignore_net;
		}

		// Set up the velocity constraint
		if (!(i_params ["equations.pressure.ignore"].IsDefined () && i_params ["equations.pressure.ignore"].as <bool> ())) {
			*div <double> (equations ["pressure"], equations ["x_velocity"], equations ["z_velocity"])
			==
			0.0;
		}
		
	TRACE ("Initialized.");
	}
	
	double boussinesq_element::calculate_timestep (int i, int j, formats::virtual_file *virtual_file) {
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
				return std::min (std::abs ((x_ptr [(i + 1) * m + j] - x_ptr [(i - 1) * m + j]) / x_vel_ptr [i * m + j]) * cfl, std::abs ((z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / z_vel_ptr [i * m + j]) * cfl);
			}
		}
		return 1.0 / 0.0;
	}
} /* pisces */
