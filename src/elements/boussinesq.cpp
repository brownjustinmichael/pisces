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

		initialize ("pressure");
		initialize ("composition");
		initialize ("temperature");
		initialize ("x_velocity");
		initialize ("z_velocity");
		initialize ("density");

		// Set up the data from the input file in params
		this->template setup_from <formats::netcdf> (i_params ["input"]);
		
		this->template setup_output_from <formats::netcdf> (i_params ["output.cart"]);

		this->template setup_output_from <formats::netcdf> (i_params ["output.trans"], transformed_horizontal);

		std::shared_ptr <io::output> stat_stream =
		this->template setup_output_from <formats::netcdf> (i_params ["output.stat"], no_variables);

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
	boussinesq_element <datatype>::boussinesq_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data <datatype> &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags, bool load_diffusion) : 
	implemented_element <datatype> (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags) {
		TRACE ("Initializing...");
		x_ptr = data ("x");
		z_ptr = data ("z");
		x_vel_ptr = data ("x_velocity", real_real);
		z_vel_ptr = data ("z_velocity", real_real);
		
		cfl = i_params ["time.cfl"].as <datatype> ();

		// Set up the temperature equation
		*split_solver <datatype> (equations ["temperature"], timestep, dirichlet (1.0), dirichlet (0.0))
		+ advec <datatype> (data ["x_velocity"], data ["z_velocity"]) 
		== 
		params ["equations.temperature.diffusion"] * diff <datatype> ();

		// // Set up the composition equation
		// *split_solver <datatype> (equations ["composition"], timestep, dirichlet (1.0), dirichlet (0.0)) 
		// + advec <datatype> (data ["x_velocity"], data ["z_velocity"]) 
		// == 
		// params ["equations.composition.diffusion"] * diff <datatype> ();

		// Set up the x_velocity equation, note the missing pressure term, which is handled in div
		*split_solver <datatype> (equations ["x_velocity"], timestep, neumann (0.0), neumann (0.0)) 
		+ advec <datatype> (data ["x_velocity"], data ["z_velocity"]) 
		== 
		params ["equations.velocity.diffusion"] * diff <datatype> ();

		// Set up the z_velocity equation, note the missing pressure term, which is handled in div
		*split_solver <datatype> (equations ["z_velocity"], timestep, dirichlet (0.0), dirichlet (0.0)) 
		+ advec <datatype> (data ["x_velocity"], data ["z_velocity"]) 
		== 
		params ["equations.z_velocity.sources.temperature"] * src <datatype> (data ["temperature"])
		// + params ["equations.z_velocity.sources.composition"] * src <datatype> (data ["composition"]) 
		+ params ["equations.velocity.diffusion"] * diff <datatype> ();
		data ["z_velocity"].component_flags |= ignore_net;

		// Set up the velocity constraint
		*div <datatype> (equations ["pressure"], equations ["x_velocity"], equations ["z_velocity"])
		==
		0.0;
		
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
				return std::min (std::abs ((x_ptr [(i + 1) * m + j] - x_ptr [(i - 1) * m + j]) / x_vel_ptr [i * m + j]), std::abs ((z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / z_vel_ptr [i * m + j])) * cfl;
			}
		}
		return 1.0 / 0.0;
	}
	
	template class boussinesq_element <double>;
} /* pisces */
