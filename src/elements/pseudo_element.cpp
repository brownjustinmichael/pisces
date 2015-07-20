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


		for (int j = 0; j < m; ++j)
		{
			data ["bg_pressure"].ptr () [j] = 1.0 + (*grids [1]) [j] * 0.1;
		}

		data ["density"] == data ["composition"] / data ["temperature"] * data ["bg_pressure"];

		// Set up the x_velocity equation, note the missing pressure term, which is handled in div
		*split_solver <datatype> (equations ["x_momentum"], timestep, neumann (0.0), neumann (0.0)) 
		+ advec <datatype> (data ["x_velocity"], data ["z_velocity"]) 
		== 
		params ["equations.velocity.diffusion"] * diff <datatype> ();
		// + horizontal_stress (data ["z_velocity"]);

		// Set up the z_velocity equation, note the missing pressure term, which is handled in div
		*split_solver <datatype> (equations ["z_momentum"], timestep, dirichlet (0.0), dirichlet (0.0)) 
		+ advec <datatype> (data ["x_velocity"], data ["z_velocity"]) 
		== 
		// - grd <datatype> (data ["bg_pressure"])
		// params ["equations.z_velocity.sources.density"] * src <datatype> (data ["density"])
		params ["equations.velocity.diffusion"] * diff <datatype> ();
		// + vertical_stress (data ["x_velocity"]);

		// // Set up the velocity constraint
		// *pdiv <datatype> (equations ["pressure"], equations ["x_momentum"], equations ["z_momentum"], data ["x_velocity"], data ["z_velocity"], data ["density"], data ["bg_pressure"])
		// ==
		// 0.0;
		
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

namespace data
{
	template <class datatype>
	pseudo_data <datatype>::pseudo_data (grids::axis *i_axis_n, grids::axis *i_axis_m, int id, int n_elements, io::parameters& i_params) : implemented_data <datatype> (i_axis_n, i_axis_m, i_params, id, i_params.get <std::string> ("dump.file"), i_params.get <std::string> ("root") + i_params.get <std::string> ("dump.directory"), i_params.get <int> ("dump.every")) {
		
		initialize ("pressure");
		initialize ("composition");
		initialize ("temperature");
		initialize ("x_velocity");
		initialize ("x_momentum");
		initialize ("z_velocity");
		initialize ("z_momentum");
		initialize ("density");

		(*this) ["x_velocity"] == (*this) ["x_momentum"] / (*this) ["density"];
		(*this) ["z_velocity"] == (*this) ["z_momentum"] / (*this) ["density"];

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
	
	template class pseudo_data <double>;
} /* data */
