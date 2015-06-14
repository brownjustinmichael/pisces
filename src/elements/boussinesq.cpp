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
#include "io/functors/slice.hpp"

#include "plans/plan.hpp"
#include "plans/advection.hpp"
#include "plans/diffusion.hpp"
#include "plans/source.hpp"
#include "plans-solvers/boundaries/implemented_boundary.hpp"
#include "plans-solvers/solvers.hpp"

namespace data
{
	template <class datatype>
	thermo_compositional_data <datatype>::thermo_compositional_data (grids::axis *i_axis_n, grids::axis *i_axis_m, int id, int n_elements, io::parameters& i_params) : implemented_data <datatype> (i_axis_n, i_axis_m, id, i_params.get <std::string> ("dump.file"), i_params.get <std::string> ("root") + i_params.get <std::string> ("dump.directory"), i_params.get <int> ("dump.every")) {
		for (YAML::const_iterator iter = i_params ["equations"].begin (); iter != i_params ["equations"].end (); ++iter) {
			initialize (iter->first.as <std::string> ());
		}

		int name = id;
		
		const formats::data_grid i_grid = formats::data_grid::two_d (n, m, 0, i_params.get <bool> ("input.full") ? n_elements * m : 0, 0, i_params.get <bool> ("input.full") ? id * m : 0);
		
		// Set up the data from the input file in params
		if (i_params.get <std::string> ("input.file") != "") {
			std::string file_format = i_params.get <std::string> ("root") + i_params.get <std::string> ("input.directory") + i_params.get <std::string> ("input.file");
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);
			io::formatted_input <formats::netcdf> input_stream (i_grid, buffer);

			(*this).setup (&input_stream);
		}
		
		// For weighted averages, calculate area
		area.resize (n * m);
		for (int i = 1; i < n; ++i) {
			for (int j = 1; j < m; ++j) {
				area [i * m + j] = ((*(this->grid_n)) [i] - (*(this->grid_n)) [i - 1]) * ((*(this->grid_m)) [j] - (*(this->grid_m)) [j - 1]);
			}
		}
		
		for (YAML::const_iterator iter = i_params ["output.files"].begin (); iter != i_params ["output.files"].end (); ++iter) {
			std::string file = iter->first.as <std::string> ();
			io::parameters::alias specs (i_params, "output.files." + file);
			
			if (specs ["output"].IsDefined ()) {
				if (!(specs ["output"].as <bool> ())) {
					continue;
				}
			} else {
				if (i_params ["output.output"].IsDefined () && !(i_params ["output.output"].as <bool> ())) {
					continue;
				}
			}
			
			std::shared_ptr <io::output> stream;
			std::string file_format = i_params.get <std::string> ("root");
			if (specs ["directory"].IsDefined ()) {
				file_format += specs ["directory"].as <std::string> ();
			} else {
				file_format += i_params ["output.directory"].as <std::string> ();
			}
			if (i_params ["output.name"].as <std::string> () != "") {
				file_format += i_params ["output.name"].as <std::string> () + "_";
			}
			file_format += file;
			
			DEBUG ("File is " << file_format << " from " << file << " " << i_params ["output.name"] << " " << specs ["directory"] << " " << i_params.get <std::string> ("root"));
			
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);
			file_format = buffer;
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), i_params ["output.count"].IsDefined () ? i_params ["output.count"].as <int> () : 0);
			
			bool full = (specs ["full"].IsDefined () && specs ["full"].as <bool> ());
			const formats::data_grid o_grid = formats::data_grid::two_d (n, m, 0, full ? n_elements * m : 0, 0, full ? id * m : 0);
			
			if (specs ["timed"].IsDefined () && specs ["timed"].as <bool> ()) {
				stream.reset (new io::timed_appender_output <datatype, formats::netcdf> (o_grid, buffer, duration, specs ["every"].IsDefined () ? specs ["every"].as <datatype> () : 1.0));
			} else {
				stream.reset (new io::appender_output <formats::netcdf> (o_grid, buffer, specs ["every"].IsDefined () ? specs ["every"].as <int> () : 1));
			}
			
			if (specs ["stat"].IsDefined () and specs ["stat"].as <bool> ()) {
				for (typename data <datatype>::iterator iter = this->begin (); iter != this->end (); ++iter) {
					// For each data variable, output z_flux, average derivative across the center, average and max
					std::string variable = *iter;
					if ((*this) ("z_velocity")) {
						stream->template append <double> ("z_flux_" + variable, std::shared_ptr <functors::functor> (new functors::average_functor <double> (n, 1, std::shared_ptr <functors::functor> (new functors::slice_functor <double> (n, m, m / 2, std::shared_ptr <functors::functor> (new functors::product_functor <double> (n, m, (*this) ("z_velocity"), (*this) (variable))))))), formats::scalar);
						stream->template append <double> ("deriv_" + variable, std::shared_ptr <functors::functor> (new functors::average_functor <double> (n, 1, std::shared_ptr <functors::functor> (new functors::slice_functor <double> (n, m, m / 2, std::shared_ptr <functors::functor> (new functors::deriv_functor <double> ((*this) (variable), n, m, &(*grid_m) [0])))))), formats::scalar);
					}
					stream->template append <double> ("avg_" + variable, std::shared_ptr <functors::functor> (new functors::weighted_average_functor <double> (n, m, &area [0], (*this) (variable))), formats::scalar);
					stream->template append <double> ("max_" + variable, std::shared_ptr <functors::functor> (new functors::max_functor <double> (n, m, (*this) (variable))), formats::scalar);
				}
			}
			
			this->setup_output (stream, ((specs ["transform"].IsDefined () and specs ["transform"].as <bool> ()) ? transformed_horizontal : 0x0) | ((specs ["stat"].IsDefined () and specs ["stat"].as <bool> ()) ? no_variables : 0x0));
			// if ((*this) ("x_velocity") && (*this) ("z_velocity")) {
			// 	normal_stream->template append <double> ("div", std::shared_ptr <functors::functor> (new functors::div_functor <double> ((*this) ("x"), (*this) ("z"), (*this) ("x_velocity"), (*this) ("z_velocity"), n, m)));
			// }
		}
		/*
			TODO Setting up the streams should be more convenient
		*/
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
		x_vel_ptr = data ("x_velocity");
		z_vel_ptr = data ("z_velocity");
		
		advection_coeff = 0.0;
		cfl = i_params.get <datatype> ("time.cfl");

		// If we aren't at an edge, add the appropriate communicating boundary

		std::shared_ptr <typename boundaries::boundary <datatype>::factory> local_boundary_0, local_boundary_n;
		std::string variable;
		for (YAML::const_iterator iter = i_params ["equations"].begin (); iter != i_params ["equations"].end (); ++iter) {
			variable = iter->first.as <std::string> ();
			io::parameters::alias terms (i_params, "equations." + variable);
			
			// If the equation is labeled as ignore, do not construct an equation
			if ((terms ["ignore"].IsDefined () && terms ["ignore"].as <bool> ()) || variable == "pressure") {
				continue;
			}
			
			// If no communicating boundary is available, we're on an edge, so use the appropriate edge case
			if (!local_boundary_0) {
				if (terms ["bottom"].IsDefined ()) {
					if (terms ["bottom.type"].as <std::string> () == "fixed_value") {
						local_boundary_0 = std::shared_ptr <typename boundaries::boundary <datatype>::factory> (new typename boundaries::fixed_boundary <datatype>::factory (terms ["bottom.value"].as <datatype> ()));
					} else if (terms ["bottom.type"].as <std::string> () == "fixed_derivative") {
						local_boundary_0 = std::shared_ptr <typename boundaries::boundary <datatype>::factory> (new typename boundaries::fixed_deriv_boundary <datatype>::factory (terms ["bottom.value"].as <datatype> ()));
					}
				}
			}
			if (!local_boundary_n) {
				if (terms ["top"].IsDefined ()) {
					if (terms ["top.type"].as <std::string> () == "fixed_value") {
						local_boundary_n = std::shared_ptr <typename boundaries::boundary <datatype>::factory> (new typename fixed_boundary <datatype>::factory (terms ["top.value"].as <datatype> ()));
					} else if (terms ["top.type"].as <std::string> () == "fixed_derivative") {
						local_boundary_n = std::shared_ptr <typename boundaries::boundary <datatype>::factory> (new typename fixed_deriv_boundary <datatype>::factory (terms ["top.value"].as <datatype> ()));
					}
				}
			}
			
			// Add the split directional solvers
			*plans::split_solver <datatype> (equations [variable],  timestep, local_boundary_0, local_boundary_n) + terms ["advection"] * advec <datatype> (ptr ("x_velocity"), ptr ("z_velocity")) == terms ["diffusion"] * diff <datatype> ();
			// equations [variable]->add_solver (typename collocation <datatype>::factory (messenger_ptr, timestep, local_boundary_0, local_boundary_n), z_solver);
			// equations [variable]->add_solver (typename fourier <datatype>::factory (timestep, local_boundary_0, local_boundary_n), x_solver);
			
			// *equations [variable] + terms ["advection"] * advec <datatype> (ptr ("x_velocity"), ptr ("z_velocity")) == terms ["diffusion"] * diff <datatype> (ptr (variable));
			
			// If any source terms are specified, construct the appropriate source plans
			if (terms ["sources"].IsDefined ()) {
				for (YAML::const_iterator source_iter = terms ["sources"].begin (); source_iter != terms ["sources"].end (); ++source_iter) {
					if (ptr (source_iter->first.as <std::string> ())) {
						*equations [variable] == source_iter->second * src (ptr (source_iter->first.as <std::string> ()));
					}
				}
			}
		}

		// *split_solver <datatype> (equations ["temperature"], timestep, dirichlet (1.0), dirichlet (0.0)) + advec <datatype> (ptr ("x_velocity"), ptr ("z_velocity")) == params ["temperature.diffusion"] * diff <datatype> ();

		// *split_solver <datatype> (equations ["composition"], timestep, dirichlet (1.0), dirichlet (0.0)) + advec <datatype> (ptr ("x_velocity"), ptr ("z_velocity")) == params ["composition.diffusion"] * diff <datatype> ();

		// *split_solver <datatype> (equations ["x_velocity"], timestep, neumann (0.0), neumann (0.0)) + advec <datatype> (ptr ("x_velocity"), ptr ("z_velocity")) == params ["velocity.diffusion"] * diff <datatype> ();

		// *split_solver <datatype> (equations ["z_velocity"], timestep, dirichlet (0.0), dirichlet (0.0)) + advec <datatype> (ptr ("x_velocity"), ptr ("z_velocity")) == params ["z_velocity.sources.temperature"] * src <datatype> (ptr ("temperature")) + params ["z_velocity.sources.composition"] * src <datatype> (ptr ("composition")) + params ["velocity.diffusion"] * diff <datatype> ();

		div <datatype> (equations ["pressure"], equations ["x_velocity"], equations ["z_velocity"]);
		
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
