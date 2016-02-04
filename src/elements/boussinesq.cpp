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
		// initialize ("density");

		// Set up the data from the input file in params
		TRACE ("Setting up from input...");
		i_params ["input.directory"] = i_params ["root"].as <std::string> () + "/" + i_params ["input.directory"].as <std::string> ();
		this->setup_from <formats::netcdf> (i_params ["input"]);
		
		i_params ["output.directory"] = i_params ["root"].as <std::string> () + "/" + i_params ["output.directory"].as <std::string> ();

		TRACE ("Setting up cartesian output...");
		std::shared_ptr <io::output> stream = this->setup_output_from <formats::netcdf> (i_params ["output.cart"]);
		// stream->append <double> ("div_u", std::shared_ptr <functors::div_functor <double>> (new functors::div_functor <double> ((*this) ["x"].ptr (), (*this) ["z"].ptr (), (*this) ["x_velocity"].ptr (), (*this) ["z_velocity"].ptr (), n, m)));

		TRACE ("Setting up transformed output...");
		stream = this->setup_output_from <formats::netcdf> (i_params ["output.trans"], real_spectral);
		// stream->append <double> ("div_u", std::shared_ptr <functors::transform_div_functor <double>> (new functors::transform_div_functor <double> ((*this) ["x"].ptr (), (*this) ["z"].ptr (), (*this) ["x_velocity"].ptr (), (*this) ["z_velocity"].ptr (), n, m)));

		TRACE ("Setting up stat output...");
		std::shared_ptr <io::output> stat_stream =
		this->setup_output_from <formats::netcdf> (i_params ["output.stat"], real_real, no_variables);

		if (!stat_stream) return;

		for (data::iterator iter = this->begin (); iter != this->end (); ++iter) {
			// For each data variable, output z_flux, average derivative across the center, average and max
			std::string variable = *iter;
			stat_stream->append <double> ("max_" + variable, this->output_max (variable), formats::scalar);
			stat_stream->append <double> ("avg_" + variable, this->output_avg (variable), formats::scalar);
			stat_stream->append <double> ("deriv_" + variable, this->output_deriv (variable), formats::scalar);
			stat_stream->append <double> ("flux_" + variable, this->output_flux (variable, "z_velocity"), formats::scalar);
		}
	}
} /* data */

namespace pisces
{
	using namespace plans;
	using namespace plans::solvers;
	using namespace boundaries;

	std::string boussinesq_element::class_name() {
		return "boussinesq";
	}

	implemented_element::registrar <boussinesq_element> boussinesq_registrar ("boussinesq");

	std::shared_ptr <element> boussinesq_element::instance (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) {
		return std::shared_ptr <element> (new boussinesq_element (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags));
	}
	
	boussinesq_element::boussinesq_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) : 
	implemented_element (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags) {
		TRACE ("Initializing...");
		x_ptr = data ("x");
		z_ptr = data ("z");
		x_vel_ptr = data ("x_velocity", real_real);
		z_vel_ptr = data ("z_velocity", real_real);
		
		cfl = i_params ["time.cfl"].as <double> ();
		allow = i_params ["time.allow"].as <double> ();
		
		// Set up the temperature equation
		if (!(i_params ["equations.temperature.ignore"].IsDefined () && i_params ["equations.temperature.ignore"].as <bool> ())) {
			*split_solver (equations ["temperature"], timestep, 
				dirichlet (i_params ["equations.temperature.bottom.value"].as <double> ()), 
				dirichlet (i_params ["equations.temperature.top.value"].as <double> ())) 
			+ params ["equations.temperature.advection"] * advec (data ["x_velocity"], data ["z_velocity"])
			// + src (data ["z_velocity"] * data ["korre_Ts"])
			== 
			params ["equations.temperature.sources.z_velocity"] * src (data ["z_velocity"])
			+ params ["equations.temperature.diffusion"] * diff ();

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
	
	/*!**********************************************************************
	 * \copydoc implemented_element::make_virtual_file
	 ************************************************************************/
	formats::virtual_file* boussinesq_element::make_virtual_file (int flags) {
		std::shared_ptr <io::output> virtual_output;
		if (flags & profile_only) {
			// If only the profile is desired, just build that
			virtual_output.reset (new io::formatted_output <formats::virtual_format> (formats::data_grid::two_d (1, m), "two_d/boussinesq/virtual_file", formats::replace_file));
			if (flags & timestep_only) {
				// If only the timestep is needed, only load z and x_velocity
				virtual_output->append <double> ("z", std::shared_ptr <functors::functor> (new functors::average_functor <double> (ptr ("z"), n, m)));
				virtual_output->append <double> ("z_velocity", std::shared_ptr <functors::functor> (new functors::root_mean_square_functor <double> (ptr ("z_velocity"), n, m)));
				virtual_output->append <double> ("temperature", std::shared_ptr <functors::functor> (new functors::root_mean_square_functor <double> (ptr ("temperature"), n, m)));
				virtual_output->append <double> ("composition", std::shared_ptr <functors::functor> (new functors::root_mean_square_functor <double> (ptr ("composition"), n, m)));
			} else {
				FATAL ("HAVEN'T GOT A TREATMENT FOR THIS YET");
				throw 0;
				data.setup_profile (virtual_output, data::no_save);
			}
		} else {
			virtual_output.reset (new io::formatted_output <formats::virtual_format> (formats::data_grid::two_d (n, m), "two_d/boussinesq/virtual_file", formats::replace_file));
			if (flags & timestep_only) {
				// If only the timestep is needed, grab the positions and velocities
				virtual_output->append <double> ("z", ptr ("z"));
				virtual_output->append <double> ("x", ptr ("x"));
				virtual_output->append <double> ("z_velocity", ptr ("z_velocity"));
				virtual_output->append <double> ("x_velocity", ptr ("x_velocity"));
				virtual_output->append <double> ("temperature", ptr ("temperature"));
				virtual_output->append <double> ("composition", ptr ("composition"));
			} else {
				// Load the whole dataset
				data.setup_output (virtual_output, data::no_save);
			}
		}
		virtual_output->to_file ();
		return &formats::virtual_files ["two_d/boussinesq/virtual_file"];
	}
	
	/*!**********************************************************************
	 * \copydoc implemented_element::make_rezoned_virtual_file
	 ************************************************************************/
	formats::virtual_file* boussinesq_element::make_rezoned_virtual_file (double *positions, formats::virtual_file *virtual_file_ptr, int flags) {
		grids::axis vertical_axis (m, positions [messenger_ptr->get_id ()], positions [messenger_ptr->get_id () + 1], messenger_ptr->get_id () == 0 ? 0 : 1, messenger_ptr->get_id () == messenger_ptr->get_np () - 1 ? 0 : 1);
		std::shared_ptr <grids::grid> vertical_grid = implemented_element::generate_grid (&vertical_axis);
		
		pisces::rezone (messenger_ptr, &*(grids [1]), &*vertical_grid, virtual_file_ptr, &formats::virtual_files ["two_d/boussinesq/new_virtual_file"], value_buffer, inter_buffer);
		
		return &formats::virtual_files ["two_d/boussinesq/new_virtual_file"];
	}
} /* pisces */
