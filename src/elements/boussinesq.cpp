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
#include "plans-solvers/solvers/collocation.hpp"
#include "plans-solvers/solvers/fourier.hpp"
#include "plans-solvers/solvers/incompressible.hpp"

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
		
		DEBUG (i_params ["output"]);
		
		for (YAML::const_iterator iter = i_params ["output.files"].begin (); iter != i_params ["output.files"].end (); ++iter) {
			std::string file = iter->first.as <std::string> ();
			io::parameters::alias specs (i_params, "output.files." + file);
			DEBUG ("Reading specs for " << file);
			
			std::shared_ptr <io::output> stream;
			std::string file_format = i_params.get <std::string> ("root");
			if (specs ["directory"].IsDefined ()) {
				file_format += specs ["directory"].as <std::string> ();
			} else {
				file_format += i_params ["output.directory"].as <std::string> ();
			}
			file_format += i_params ["output.name"].as <std::string> () + "_" + file;
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);
			snprintf (buffer, file_format.size () * 2, buffer, i_params ["output.count"].IsDefined () ? i_params ["output.count"].as <int> () : 0);
			
			bool full = (specs ["full"].IsDefined () && specs ["full"].as <bool> ());
			const formats::data_grid o_grid = formats::data_grid::two_d (n, m, 0, full ? n_elements * m : 0, 0, full ? id * m : 0);
			
			if (specs ["timed"].IsDefined () && specs ["timed"].as <bool> ()) {
				stream.reset (new io::timed_appender_output <datatype, formats::netcdf> (o_grid, buffer, duration, specs ["every"].IsDefined () ? specs ["every"].as <datatype> () : 1.0));
			} else {
				stream.reset (new io::appender_output <formats::netcdf> (o_grid, buffer, specs ["every"].IsDefined () ? specs ["every"].as <int> () : 1));
			}
			
			if (specs ["stat"].IsDefined () and specs ["stat"].as <bool> ()) {
				this->setup_stat (stream);
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
			if (specs ["transform"].IsDefined () and specs ["transform"].as <bool> ()) {
				this->setup_output (stream, transformed_horizontal);
			} else {
				this->setup_output (stream);
			}
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
		std::shared_ptr <boundaries::boundary <datatype>> boundary_0, boundary_n;
		if (messenger_ptr->get_id () > 0) {
			boundary_0 = std::shared_ptr <boundaries::boundary <datatype>> (new communicating_boundary <datatype> (messenger_ptr, grids [0]->get_ld (), m, grids [1]->get_excess_0 (), false));
		}
		if (messenger_ptr->get_id () + 1 < messenger_ptr->get_np ()) {
			boundary_n = std::shared_ptr <boundaries::boundary <datatype>> (new communicating_boundary <datatype> (messenger_ptr, grids [0]->get_ld (), m, grids [1]->get_excess_n (), true));
		}

		std::shared_ptr <boundaries::boundary <datatype>> local_boundary_0, local_boundary_n;
		std::string variable;
		for (YAML::const_iterator iter = i_params ["equations"].begin (); iter != i_params ["equations"].end (); ++iter) {
			variable = iter->first.as <std::string> ();
			io::parameters::alias terms (i_params, "equations." + variable);
			
			// If the equation is labeled as ignore, do not construct an equation
			if ((terms ["ignore"].IsDefined () && terms ["ignore"].as <bool> ()) || variable == "pressure") {
				continue;
			}
			
			// Use any available communicating boundaries
			local_boundary_0 = boundary_0;
			local_boundary_n = boundary_n;
			
			// If no communicating boundary is available, we're on an edge, so use the appropriate edge case
			if (!local_boundary_0) {
				if (terms ["bottom"].IsDefined ()) {
					if (terms ["bottom.type"].as <std::string> () == "fixed_value") {
						local_boundary_0 = std::shared_ptr <boundaries::boundary <datatype>> (new fixed_boundary <datatype> (&*grids [0], &*grids [1], terms ["bottom.value"].as <datatype> (), false));
					} else if (terms ["bottom.type"].as <std::string> () == "fixed_derivative") {
						local_boundary_0 = std::shared_ptr <boundaries::boundary <datatype>> (new fixed_deriv_boundary <datatype> (&*grids [0], &*grids [1], terms ["bottom.value"].as <datatype> (), false));
					}
				}
			}
			if (!local_boundary_n) {
				if (terms ["top"].IsDefined ()) {
					if (terms ["top.type"].as <std::string> () == "fixed_value") {
						local_boundary_n = std::shared_ptr <boundaries::boundary <datatype>> (new fixed_boundary <datatype> (&*grids [0], &*grids [1], terms ["top.value"].as <datatype> (), true));
					} else if (terms ["top.type"].as <std::string> () == "fixed_derivative") {
						local_boundary_n = std::shared_ptr <boundaries::boundary <datatype>> (new fixed_deriv_boundary <datatype> (&*grids [0], &*grids [1], terms ["top.value"].as <datatype> (), true));
					}
				}
			}
			
			// Add the split directional solvers
			equations [variable]->add_solver (typename collocation <datatype>::factory (messenger_ptr, timestep, local_boundary_0, local_boundary_n), z_solver);
			equations [variable]->add_solver (typename fourier <datatype>::factory (timestep, local_boundary_0, local_boundary_n), x_solver);
			
			if (load_diffusion) {
				// If a diffusion value is specified, construct the diffusion plans
				equations [variable]->add_plan (typename diffusion::vertical <datatype>::factory (terms ["diffusion"], i_params.get <datatype> ("time.alpha")), pre_plan);
				equations [variable]->add_plan (typename diffusion::horizontal <datatype>::factory (terms ["diffusion"], i_params.get <datatype> ("time.alpha")), mid_plan);
			}
			
			// If an advection value is specified, construct the advection plan
			equations [variable]->add_plan (typename advection::uniform <datatype>::factory (terms ["advection"], ptr ("x_velocity"), ptr ("z_velocity")), post_plan);
			if (terms ["advection"].IsDefined ()) advection_coeff = std::max (advection_coeff, terms ["advection"].as <datatype> ());
			
			// If any source terms are specified, construct the appropriate source plans
			if (terms ["sources"].IsDefined ()) {
				for (YAML::const_iterator source_iter = terms ["sources"].begin (); source_iter != terms ["sources"].end (); ++source_iter) {
					if (ptr (source_iter->first.as <std::string> ())) {
						equations [variable]->add_plan (typename source::uniform <datatype>::factory (source_iter->second.as <datatype> (), ptr (source_iter->first.as <std::string> ())), mid_plan);
					}
				}
			}
		}
		
		// Since this is a Boussinesq problem, also include the pressure term
		if (!i_params ["equations.pressure.ignore"].as <bool> ()) {
			if (!ptr (i_params ["equations.pressure.x_velocity"].as <std::string> ()) || !ptr (i_params ["equations.pressure.z_velocity"].as <std::string> ())) {
				FATAL ("Uninitialized velocity specified in pressure solve");
			}
			equations ["pressure"]->add_solver (typename incompressible <datatype>::factory (messenger_ptr, boundary_0, boundary_n, *equations ["x_velocity"], *equations ["z_velocity"]));
			equations ["pressure"]->get_solver (x_solver)->add_dependency ("z_velocity");
			equations ["pressure"]->get_solver (x_solver)->add_dependency ("x_velocity");
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
				return std::abs ((virtual_file->index <datatype> ("z", i, j + 1) - virtual_file->index <datatype> ("z", i, j - 1)) / virtual_file->index <datatype> ("z_velocity", i, j) / advection_coeff) * cfl;
			} else {
				return std::min (std::abs ((virtual_file->index <datatype> ("x", i + 1, j) - virtual_file->index <datatype> ("x", i - 1, j)) / virtual_file->index <datatype> ("x_velocity", i, j) / advection_coeff), std::abs ((virtual_file->index <datatype> ("z", i, j + 1) - virtual_file->index <datatype> ("z", i, j - 1)) / virtual_file->index <datatype> ("z_velocity", i, j) / advection_coeff)) * cfl;
			}
		} else {
			if (j == 0 || j == m - 1) {
				return 1.0 / 0.0;
			}
			if ((i == 0 || i == n - 1)) {
				return std::abs ((z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / z_vel_ptr [i * m + j] / advection_coeff) * cfl;
			} else {
				return std::min (std::abs ((x_ptr [(i + 1) * m + j] - x_ptr [(i - 1) * m + j]) / x_vel_ptr [i * m + j] / advection_coeff), std::abs ((z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / z_vel_ptr [i * m + j] / advection_coeff)) * cfl;
			}
		}
		return 1.0 / 0.0;
	}
	
	template class boussinesq_element <double>;
} /* pisces */
