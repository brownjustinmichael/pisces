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
		initialize (x_velocity, "u");
		initialize (z_velocity, "w");
		initialize (temp, "T");
		initialize (composition, "S");
		initialize (pressure, "P");
		
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
			normal_stream->template append <double> ("div", new io::functors::div_functor <double> ((*this) (x_position), (*this) (z_position), (*this) (x_velocity), (*this) (z_velocity), n, m));
		}

		std::shared_ptr <io::output> transform_stream;
		if (i_params.get <std::string> ("output.transform.file") != "") {
			std::string file_format = i_params.get <std::string> ("root") + i_params.get <std::string> ("output.directory") + i_params.get <std::string> ("output.transform.file");
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);

			transform_stream.reset (new io::appender_output <io::formats::two_d::netcdf> (o_grid, buffer, i_params.get <int> ("output.every")));
			this->setup_output (transform_stream, transformed_horizontal);
			transform_stream->template append <double> ("div", new io::functors::transform_div_functor <double> ((*this) (x_position), (*this) (z_position), (*this) (x_velocity), (*this) (z_velocity), n, m));
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
			stat_stream->template append <double> ("wT", new io::functors::weighted_average_functor <double> (n, m, &area [0], new io::functors::product_functor <double> (n, m, (*this) (z_velocity), (*this) (temperature))), io::scalar);
			stat_stream->template append <double> ("wS", new io::functors::weighted_average_functor <double> (n, m, &area [0], new io::functors::product_functor <double> (n, m, (*this) (z_velocity), (*this) (composition))), io::scalar);
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
		x_ptr = data (x_position);
		z_ptr = data (z_position);
		x_vel_ptr = data (x_velocity);
		z_vel_ptr = data (z_velocity);
		
		advection_coeff = i_params.get <datatype> ("temperature.advection");
		advection_coeff = std::max (advection_coeff, i_params.get <datatype> ("velocity.advection"));
		advection_coeff = std::max (advection_coeff, i_params.get <datatype> ("composition.advection"));
		cfl = i_params.get <datatype> ("time.cfl");

		std::shared_ptr <plans::boundary <datatype>> boundary_0, boundary_n, deriv_boundary_0, deriv_boundary_n, thermal_boundary_0, thermal_boundary_n, compositional_boundary_0, compositional_boundary_n;
		if (messenger_ptr->get_id () > 0) {
			boundary_0 = std::shared_ptr <plans::boundary <datatype>> (new communicating_boundary <datatype> (messenger_ptr, grids [0]->get_ld (), m, grids [1]->get_excess_0 (), &((*grids [1]) [0]), 0, false));
			deriv_boundary_0 = boundary_0;
			thermal_boundary_0 = boundary_0;
			compositional_boundary_0 = boundary_0;
		} else {
			boundary_0 = std::shared_ptr <plans::boundary <datatype>> (new fixed_boundary <datatype> (&*grids [0], &*grids [1], 0.0, false));
			deriv_boundary_0 = std::shared_ptr <plans::boundary <datatype>> (new fixed_deriv_boundary <datatype> (&*grids [0], &*grids [1], 0.0, false));
			// thermal_boundary_0 = deriv_boundary_0;
			thermal_boundary_0 = boundary_0;
			// compositional_boundary_0 = deriv_boundary_0;
			compositional_boundary_0 = boundary_0;
		}
		if (messenger_ptr->get_id () + 1 < messenger_ptr->get_np ()) {
			boundary_n = std::shared_ptr <plans::boundary <datatype>> (new communicating_boundary <datatype> (messenger_ptr, grids [0]->get_ld (), m, grids [1]->get_excess_n (), &((*grids [1]) [0]), m - grids [1]->get_excess_n (), true));
			deriv_boundary_n = boundary_n;
			thermal_boundary_n = boundary_n;
			compositional_boundary_n = boundary_n;
		} else {
			boundary_n = std::shared_ptr <plans::boundary <datatype>> (new fixed_boundary <datatype> (&*grids [0], &*grids [1], 0.0, true));
			deriv_boundary_n = std::shared_ptr <plans::boundary <datatype>> (new fixed_deriv_boundary <datatype> (&*grids [0], &*grids [1], 0.0, true));
			// thermal_boundary_n = deriv_boundary_n;
			thermal_boundary_n = boundary_n;
			// compositional_boundary_n = deriv_boundary_n;
			compositional_boundary_n = boundary_n;
		}

		/*
			TODO Figure out how to more conveniently determine whether an edge effect is needed.
		*/
		
		// diffusion.resize (m);
		// for (int j = 0; j < m; ++j) {
		// 	diffusion [j] = (*grids [1]) [j] > 0.0 ? i_params.get <datatype> ("temperature.diffusion") : i_params.get <datatype> ("temperature.bg_diffusion");
		// }
		
		// Solve temperature
		solvers [temp]->add_solver (typename collocation_solver <datatype>::factory (messenger_ptr, timestep, thermal_boundary_0, thermal_boundary_n), z_solver);
		solvers [temp]->add_solver (typename fourier_solver <datatype>::factory (timestep, thermal_boundary_0, thermal_boundary_n), x_solver);

		solvers [temp]->add_plan (typename vertical_diffusion <datatype>::factory (i_params.get <datatype> ("temperature.diffusion"), i_params.get <datatype> ("time.alpha")), pre_plan);
		solvers [temp]->add_plan (typename horizontal_diffusion <datatype>::factory (i_params.get <datatype> ("temperature.diffusion"), i_params.get <datatype> ("time.alpha")), mid_plan);
		solvers [temp]->add_plan (typename advection <datatype>::factory (i_params.get <datatype> ("temperature.advection"), ptr (x_vel), ptr (z_vel)), post_plan);
		solvers [temp]->add_plan (typename source <datatype>::factory (-i_params.get <datatype> ("temperature.stratification"), ptr (z_velocity)), mid_plan);

		// Solve composition
		solvers [composition]->add_solver (typename collocation_solver <datatype>::factory (messenger_ptr, timestep, compositional_boundary_0, compositional_boundary_n), z_solver);
		solvers [composition]->add_solver (typename fourier_solver <datatype>::factory (timestep, compositional_boundary_0, compositional_boundary_n), x_solver);

		solvers [composition]->add_plan (typename vertical_diffusion <datatype>::factory (i_params.get <datatype> ("composition.diffusion"), i_params.get <datatype> ("time.alpha")), pre_plan);
		solvers [composition]->add_plan (typename horizontal_diffusion <datatype>::factory (i_params.get <datatype> ("composition.diffusion"), i_params.get <datatype> ("time.alpha")), mid_plan);
		solvers [composition]->add_plan (typename advection <datatype>::factory (i_params.get <datatype> ("composition.advection"), ptr (x_vel), ptr (z_vel)), post_plan);
		solvers [composition]->add_plan (typename source <datatype>::factory (-i_params.get <datatype> ("composition.stratification"), ptr (z_velocity)), mid_plan);
	
		// Solve velocity
		solvers [x_velocity]->add_solver (typename collocation_solver <datatype>::factory (messenger_ptr, timestep, deriv_boundary_0, deriv_boundary_n), z_solver);
		solvers [x_velocity]->add_solver (typename fourier_solver <datatype>::factory (timestep, deriv_boundary_0, deriv_boundary_n), x_solver);

		solvers [x_velocity]->add_plan (typename vertical_diffusion <datatype>::factory (i_params.get <datatype> ("velocity.diffusion"), i_params.get <datatype> ("time.alpha")), pre_plan);
		solvers [x_velocity]->add_plan (typename horizontal_diffusion <datatype>::factory (i_params.get <datatype> ("velocity.diffusion"), i_params.get <datatype> ("time.alpha")), mid_plan);
		solvers [x_velocity]->add_plan (typename advection <datatype>::factory (i_params.get <datatype> ("velocity.advection"), ptr (x_velocity), ptr (z_velocity)), post_plan);

		solvers [z_velocity]->add_solver (typename collocation_solver <datatype>::factory (messenger_ptr, timestep, boundary_0, boundary_n), z_solver);
		solvers [z_velocity]->add_solver (typename fourier_solver <datatype>::factory (timestep, boundary_0, boundary_n), x_solver);

		solvers [z_velocity]->add_plan (typename vertical_diffusion <datatype>::factory (i_params.get <datatype> ("velocity.diffusion"), i_params.get <datatype> ("time.alpha")), pre_plan);
		solvers [z_velocity]->add_plan (typename horizontal_diffusion <datatype>::factory (i_params.get <datatype> ("velocity.diffusion"), i_params.get <datatype> ("time.alpha")), mid_plan);
		solvers [z_velocity]->add_plan (typename advection <datatype>::factory (i_params.get <datatype> ("velocity.advection"), ptr (x_velocity), ptr (z_velocity)), post_plan);
		solvers [z_velocity]->add_plan (typename source <datatype>::factory (i_params.get <datatype> ("velocity.buoyancy.temperature"), ptr (temp)), mid_plan);
		solvers [z_velocity]->add_plan (typename source <datatype>::factory (i_params.get <datatype> ("velocity.buoyancy.composition"), ptr (composition)), mid_plan);

		solvers [pressure]->add_solver (typename incompressible_corrector <datatype>::factory (messenger_ptr, boundary_0, boundary_n, *solvers [x_velocity], *solvers [z_velocity]));
		solvers [pressure]->get_solver (x_solver)->add_dependency (z_velocity);
		solvers [pressure]->get_solver (x_solver)->add_dependency (x_velocity);

		TRACE ("Initialized.");
	}
	
	template <class datatype>
	datatype boussinesq_element <datatype>::calculate_timestep (int i, int j, io::formats::virtual_file *virtual_file) {
		if (virtual_file) {
			if (j == 0 || j == virtual_file->dims ["z"] [1] - 1) {
				return 1.0 / 0.0;
			}
			if (i == 0 || i == virtual_file->dims ["z"] [0] - 1) {
				return std::abs ((virtual_file->index <datatype> ("z", i, j + 1) - virtual_file->index <datatype> ("z", i, j - 1)) / virtual_file->index <datatype> ("w", i, j) / advection_coeff) * cfl;
			} else {
				return std::min (std::abs ((virtual_file->index <datatype> ("x", i + 1, j) - virtual_file->index <datatype> ("x", i - 1, j)) / virtual_file->index <datatype> ("u", i, j) / advection_coeff), std::abs ((virtual_file->index <datatype> ("z", i, j + 1) - virtual_file->index <datatype> ("z", i, j - 1)) / virtual_file->index <datatype> ("w", i, j) / advection_coeff)) * cfl;
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