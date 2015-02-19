/*!**********************************************************************
 * \file boussinesq_two_d.cpp
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
#include "plans-solvers/implemented_boundary.hpp"
#include "plans-solvers/solvers/collocation.hpp"
#include "plans-solvers/solvers/fourier.hpp"
#include "plans-solvers/solvers/incompressible.hpp"

namespace data
{
	template <class datatype>
	thermo_compositional_data <datatype>::thermo_compositional_data (plans::axis *i_axis_n, plans::axis *i_axis_m, int id, int n_elements, io::parameters& i_params) : implemented_data <datatype> (i_axis_n, i_axis_m, id, i_params.get <std::string> ("dump.file"), i_params.get <std::string> ("root") + i_params.get <std::string> ("dump.directory"), i_params.get <int> ("dump.every")) {
		for (YAML::const_iterator iter = i_params ["equations"].begin (); iter != i_params ["equations"].end (); ++iter) {
			if (!(iter->second ["ignore"].IsDefined () && iter->second ["ignore"].as <bool> ())) initialize (iter->first.as <std::string> ());
		}
		initialize ("pressure");
		
		int name = id;
		
		const io::data_grid i_grid = io::data_grid::two_d (n, m, 0, i_params.get <bool> ("input.full") ? n_elements * m : 0, 0, i_params.get <bool> ("input.full") ? id * m : 0);

		if (i_params.get <std::string> ("input.file") != "") {
			std::string file_format = i_params.get <std::string> ("root") + i_params.get <std::string> ("input.directory") + i_params.get <std::string> ("input.file");
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);
			io::formatted_input <io::formats::two_d::netcdf> input_stream (i_grid, buffer);

			(*this).setup (&input_stream);
		}

		const io::data_grid o_grid = io::data_grid::two_d (n, m, 0, i_params.get <bool> ("output.full") ? n_elements * m : 0, 0, i_params.get <bool> ("output.full") ? id * m : 0);
		
		// Set up output
		std::shared_ptr <io::output> normal_stream;
		if (i_params.get <std::string> ("output.file") != "") {
			std::string file_format = i_params.get <std::string> ("root") + i_params.get <std::string> ("output.directory") + i_params.get <std::string> ("output.file");
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);

			normal_stream.reset (new io::appender_output <io::formats::two_d::netcdf> (o_grid, buffer, i_params.get <int> ("output.every")));
			this->setup_output (normal_stream);
			normal_stream->template append <double> ("div", std::shared_ptr <io::functors::functor> (new io::functors::div_functor <double> ((*this) ("x"), (*this) ("z"), (*this) ("x_velocity"), (*this) ("z_velocity"), n, m)));
		}

		std::shared_ptr <io::output> transform_stream;
		if (i_params.get <std::string> ("output.transform.file") != "") {
			std::string file_format = i_params.get <std::string> ("root") + i_params.get <std::string> ("output.directory") + i_params.get <std::string> ("output.transform.file");
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);

			transform_stream.reset (new io::appender_output <io::formats::two_d::netcdf> (o_grid, buffer, i_params.get <int> ("output.every")));
			this->setup_output (transform_stream, transformed_horizontal);
			transform_stream->template append <double> ("div", std::shared_ptr <io::functors::functor> (new io::functors::transform_div_functor <double> ((*this) ("x"), (*this) ("z"), (*this) ("x_velocity"), (*this) ("z_velocity"), n, m)));
		}

		std::shared_ptr <io::output> stat_stream;
		if (i_params.get <std::string> ("output.stat.file") != "") {
			std::string file_format = i_params.get <std::string> ("root") + i_params.get <std::string> ("output.directory") + i_params.get <std::string> ("output.stat.file");
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);

			area.resize (n * m);
			for (int i = 1; i < n; ++i) {
				for (int j = 1; j < m; ++j) {
					area [i * m + j] = ((*(this->grid_n)) [i] - (*(this->grid_n)) [i - 1]) * ((*(this->grid_m)) [j] - (*(this->grid_m)) [j - 1]);
				}
			}

			stat_stream.reset (new io::appender_output <io::formats::ascii> (io::data_grid::two_d (n, m), buffer, i_params.get <int> ("output.stat.every")));
			this->setup_stat (stat_stream);
			stat_stream->template append <double> ("wT", std::shared_ptr <io::functors::functor> (new io::functors::weighted_average_functor <double> (n, m, &area [0], std::shared_ptr <io::functors::functor> (new io::functors::product_functor <double> (n, m, (*this) ("z_velocity"), (*this) ("temperature"))))), io::scalar);
			stat_stream->template append <double> ("wC", std::shared_ptr <io::functors::functor> (new io::functors::weighted_average_functor <double> (n, m, &area [0], std::shared_ptr <io::functors::functor> (new io::functors::product_functor <double> (n, m, (*this) ("z_velocity"), (*this) ("composition"))))), io::scalar);
			stat_stream->template append <double> ("Tavg", std::shared_ptr <io::functors::functor> (new io::functors::weighted_average_functor <double> (n, m, &area [0], (*this) ("temperature"))), io::scalar);
			stat_stream->template append <double> ("Cavg", std::shared_ptr <io::functors::functor> (new io::functors::weighted_average_functor <double> (n, m, &area [0], (*this) ("composition"))), io::scalar);
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
	
	template <class datatype>
	boussinesq_element <datatype>::boussinesq_element (plans::axis i_axis_n, plans::axis i_axis_m, int i_name, io::parameters& i_params, data::data <datatype> &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) : 
	implemented_element <datatype> (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags) {
		TRACE ("Initializing...");
		x_ptr = data ("x");
		z_ptr = data ("z");
		x_vel_ptr = data ("x_velocity");
		z_vel_ptr = data ("z_velocity");
		
		advection_coeff = 0.0;
		cfl = i_params.get <datatype> ("time.cfl");

		// If we aren't at an edge, add the appropriate communicating boundary
		std::shared_ptr <plans::boundary <datatype>> boundary_0, boundary_n;
		if (messenger_ptr->get_id () > 0) {
			boundary_0 = std::shared_ptr <plans::boundary <datatype>> (new communicating_boundary <datatype> (messenger_ptr, grids [0]->get_ld (), m, grids [1]->get_excess_0 (), &((*grids [1]) [0]), 0, false));
		}
		if (messenger_ptr->get_id () + 1 < messenger_ptr->get_np ()) {
			boundary_n = std::shared_ptr <plans::boundary <datatype>> (new communicating_boundary <datatype> (messenger_ptr, grids [0]->get_ld (), m, grids [1]->get_excess_n (), &((*grids [1]) [0]), m - grids [1]->get_excess_n (), true));
		}

		std::shared_ptr <plans::boundary <datatype>> local_boundary_0, local_boundary_n;
		std::string variable;
		for (YAML::const_iterator iter = i_params ["equations"].begin (); iter != i_params ["equations"].end (); ++iter) {
			variable = iter->first.as <std::string> ();
			io::parameters::alias terms (i_params, "equations." + variable);
			
			DEBUG (variable);
			
			// If the equation is labeled as ignore, do not construct an equation
			if (terms ["ignore"].IsDefined () && terms ["ignore"].as <bool> ()) {
				continue;
			}
			
			// Use any available communicating boundaries
			local_boundary_0 = boundary_0;
			local_boundary_n = boundary_n;
			
			// If no communicating boundary is available, we're on an edge, so use the appropriate edge case
			if (!local_boundary_0) {
				DEBUG ("Checking bottom");
				if (terms ["bottom"].IsDefined ()) {
					if (terms ["bottom.type"].as <std::string> () == "fixed_value") {
						DEBUG ("Fixed");
						local_boundary_0 = std::shared_ptr <plans::boundary <datatype>> (new fixed_boundary <datatype> (&*grids [0], &*grids [1], terms ["bottom.value"].as <datatype> (), false));
					} else if (terms ["bottom.type"].as <std::string> () == "fixed_derivative") {
						DEBUG ("Fixed deriv");
						local_boundary_0 = std::shared_ptr <plans::boundary <datatype>> (new fixed_deriv_boundary <datatype> (&*grids [0], &*grids [1], terms ["bottom.value"].as <datatype> (), false));
					}
				}
			}
			if (!local_boundary_n) {
				DEBUG ("Checking top");
				
				if (terms ["top"].IsDefined ()) {
					if (terms ["top.type"].as <std::string> () == "fixed_value") {
						DEBUG ("Fixed");
						local_boundary_0 = std::shared_ptr <plans::boundary <datatype>> (new fixed_boundary <datatype> (&*grids [0], &*grids [1], terms ["top.value"].as <datatype> (), false));
					} else if (terms ["top.type"].as <std::string> () == "fixed_derivative") {
						DEBUG ("Fixed deriv");
						local_boundary_0 = std::shared_ptr <plans::boundary <datatype>> (new fixed_deriv_boundary <datatype> (&*grids [0], &*grids [1], terms ["top.value"].as <datatype> (), false));
					}
				}
			}
			
			// Add the split directional solvers
			solvers [variable]->add_solver (typename collocation_solver <datatype>::factory (messenger_ptr, timestep, boundary_0, boundary_n), z_solver);
			solvers [variable]->add_solver (typename fourier_solver <datatype>::factory (timestep, boundary_0, boundary_n), x_solver);

			// If a diffusion value is specified, construct the diffusion plans
			solvers [variable]->add_plan (typename vertical_diffusion <datatype>::factory (terms ["diffusion"], i_params.get <datatype> ("time.alpha")), pre_plan);
			solvers [variable]->add_plan (typename horizontal_diffusion <datatype>::factory (terms ["diffusion"], i_params.get <datatype> ("time.alpha")), mid_plan);
			
			// If an advection value is specified, construct the advection plan
			solvers [variable]->add_plan (typename advection <datatype>::factory (terms ["advection"], ptr ("x_velocity"), ptr ("z_velocity")), post_plan);
			if (terms ["advection"].IsDefined ()) advection_coeff = std::max (advection_coeff, terms ["advection"].as <datatype> ());
			
			// If any source terms are specified, construct the appropriate source plans
			if (terms ["sources"].IsDefined ()) {
				for (YAML::const_iterator source_iter = terms ["sources"].begin (); source_iter != terms ["sources"].end (); ++source_iter) {
					solvers [variable]->add_plan (typename source <datatype>::factory (source_iter->second.as <datatype> (), ptr (source_iter->first.as <std::string> ())), mid_plan);
				}
			}
		}
		
		// Since this is a Boussinesq problem, also include the pressure term
		solvers ["pressure"]->add_solver (typename incompressible_corrector <datatype>::factory (messenger_ptr, boundary_0, boundary_n, *solvers ["x_velocity"], *solvers ["z_velocity"]));
		solvers ["pressure"]->get_solver (x_solver)->add_dependency ("z_velocity");
		solvers ["pressure"]->get_solver (x_solver)->add_dependency ("x_velocity");

		TRACE ("Initialized.");
	}
	
	template <class datatype>
	datatype boussinesq_element <datatype>::calculate_timestep (int i, int j, io::formats::virtual_file *virtual_file) {
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
			if (i == 0 || i == n - 1) {
				return std::abs ((z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / z_vel_ptr [i * m + j] / advection_coeff) * cfl;
			} else {
				return std::min (std::abs ((x_ptr [(i + 1) * m + j] - x_ptr [(i - 1) * m + j]) / x_vel_ptr [i * m + j] / advection_coeff), std::abs ((z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / z_vel_ptr [i * m + j] / advection_coeff)) * cfl;
			}
		}
	}
	
	template class boussinesq_element <double>;
} /* pisces */
