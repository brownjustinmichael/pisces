/*!**********************************************************************
 * \file boussinesq_two_d.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "boussinesq_two_d.hpp"
#include "advection_two_d.hpp"
#include "diffusion_two_d.hpp"
#include "transform_two_d.hpp"
#include "source_two_d.hpp"
#include "solver_two_d.hpp"
#include <sstream>
#include <iomanip>
#include <stdlib.h>

namespace two_d
{
	namespace fourier
	{
		namespace cosine
		{
			template <class datatype>
			boussinesq_element <datatype>::boussinesq_element (bases::axis *i_axis_n, bases::axis *i_axis_m, int i_name, io::parameters& i_params, bases::messenger* i_messenger_ptr, int i_element_flags) : 
			element <datatype> (i_axis_n, i_axis_m, i_name, i_params, i_messenger_ptr, i_element_flags) {
		
				assert (n > 0);
				assert (m > 0);
		
				TRACE ("Initializing...");
				x_ptr = ptr (x_position);
				z_ptr = ptr (z_position);
				initialize (temp);
				x_vel_ptr = initialize (x_velocity);
				z_vel_ptr = initialize (z_velocity);
				initialize (pressure);
				
				// Set up output
				{
					std::string file_format = "../output/" + i_params.get <std::string> ("output.file");
					char buffer [file_format.size () * 2];
					snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);

					normal_stream.reset (new io::incremental (new io::two_d::netcdf (n, m), buffer, i_params.get <int> ("output.every")));
					normal_stream->template append <int> ("i", &cell_n [0]);
					normal_stream->template append <int> ("j", &cell_m [0]);
					normal_stream->template append <datatype> ("x", ptr (x_position));
					normal_stream->template append <datatype> ("z", ptr (z_position));
					normal_stream->template append <datatype> ("u", ptr (x_velocity));
					normal_stream->template append <datatype> ("w", ptr (z_velocity));
					normal_stream->template append <datatype> ("T", ptr (temp));
					normal_stream->template append <datatype> ("P", ptr (pressure));
					normal_stream->template append_scalar <datatype> ("t", &duration);
				}
				
				{
					std::string file_format = "../output/" + i_params.get <std::string> ("output.transform_file");
					char buffer [file_format.size () * 2];
					snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);

					transform_stream.reset (new io::incremental (new io::two_d::netcdf (n, m), buffer, i_params.get <int> ("output.every")));
					transform_stream->template append <int> ("i", &cell_n [0]);
					transform_stream->template append <int> ("j", &cell_m [0]);
					transform_stream->template append <datatype> ("x", ptr (x_position));
					transform_stream->template append <datatype> ("z", ptr (z_position));
					transform_stream->template append <datatype> ("u", ptr (x_velocity));
					transform_stream->template append <datatype> ("w", ptr (z_velocity));
					transform_stream->template append <datatype> ("T", ptr (temp));
					transform_stream->template append <datatype> ("P", ptr (pressure));
					transform_stream->template append_scalar <datatype> ("t", &duration);
				}

				advection_coeff = std::max (i_params.get <datatype> ("temperature.advection"), i_params.get <datatype> ("velocity.advection"));
				cfl = i_params.get <datatype> ("time.cfl");

				// Solve velocity
				element <datatype>::add_solver (x_velocity, new divergence_solver <datatype> (*grids [0], *grids [1], ptr (x_velocity), ptr (z_velocity), &element_flags [state], &element_flags [x_velocity]));
				
				element <datatype>::add_solver (z_velocity, new solver <datatype> (*grids [0], *grids [1], messenger_ptr, timestep, alpha_0, alpha_n, ptr (z_velocity), &element_flags [state], &element_flags [z_velocity]));
				
				solvers [z_velocity]->add_pre_plan (new vertical_diffusion <datatype> (*solvers [z_velocity], i_params.get <datatype> ("velocity.diffusion"), i_params.get <datatype> ("time.alpha")));
				solvers [z_velocity]->add_mid_plan (new horizontal_diffusion <datatype> (*solvers [z_velocity], i_params.get <datatype> ("velocity.diffusion"), i_params.get <datatype> ("time.alpha")));
				solvers [z_velocity]->add_post_plan (new advection <datatype> (*solvers [z_velocity], i_params.get <datatype> ("velocity.advection"), ptr (x_velocity), ptr (z_velocity)));
				solvers [z_velocity]->add_mid_plan (new source <datatype> (*solvers [z_velocity], 1.0, ptr (temp)));
				
				std::shared_ptr <bases::solver <datatype>> pressure_solver = std::shared_ptr <bases::solver <datatype>> (new laplace_solver <datatype> (*grids [0], *grids [1], messenger_ptr, ptr (pressure), &element_flags [state], &element_flags [pressure]));
				solvers [z_velocity]->add_pre_solve_plan (pressure_solver);
				solvers [z_velocity]->add_mid_plan (new z_derivative_source <datatype> (*pressure_solver, 1.0, ptr (temp)));
				solvers [z_velocity]->add_post_plan (new mixed_derivative_source <datatype> (*pressure_solver, -2.0, ptr (z_vel), ptr (x_vel)));
				solvers [z_velocity]->add_post_plan (new square_z_derivative_source <datatype> (*pressure_solver, -1.0, ptr (z_vel)));
				solvers [z_velocity]->add_post_plan (new square_x_derivative_source <datatype> (*pressure_solver, -1.0, ptr (x_vel)));
				solvers [z_velocity]->add_pre_solve_plan (new z_derivative_source <datatype> (*solvers [z_velocity], -1.0, ptr (pressure)));
				
				// Solve temperature
				element <datatype>::add_solver (temp, new solver <datatype> (*grids [0], *grids [1], messenger_ptr, timestep, alpha_0, alpha_n, ptr (temp), &element_flags [state], &element_flags [temp]));
						
				solvers [temp]->add_pre_plan (new vertical_diffusion <datatype> (*solvers [temp], i_params.get <datatype> ("temperature.diffusion"), i_params.get <datatype> ("time.alpha")));
				solvers [temp]->add_mid_plan (new horizontal_diffusion <datatype> (*solvers [temp], i_params.get <datatype> ("temperature.diffusion"), i_params.get <datatype> ("time.alpha")));
				solvers [temp]->add_post_plan (new advection <datatype> (*solvers [temp], i_params.get <datatype> ("temperature.advection"), ptr (x_vel), ptr (z_vel)));
				solvers [temp]->add_mid_plan (new source <datatype> (*solvers [temp], -i_params.get <datatype> ("temperature.advection"), ptr (z_velocity)));
				
				// Solve density

				TRACE ("Initialized.");
			}
			
			template <class datatype>
			datatype boussinesq_element <datatype>::calculate_timestep (int i, int j) {
				return std::min (std::abs ((x_ptr [(i + 1) * m + j] - x_ptr [(i - 1) * m + j]) / x_vel_ptr [i * m + j] / advection_coeff), std::abs ((z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / z_vel_ptr [i * m + j] / advection_coeff)) * cfl;
			}
			
			template <class datatype>
			void boussinesq_element <datatype>::setup (io::input *input_stream) {
				// Set up input
				input_stream->template append <datatype> ("T", ptr (temp));
				input_stream->template append_scalar <datatype> ("t", &duration);
				int mode;
				input_stream->template append_scalar <int> ("mode", &mode);
				input_stream->from_file ();
				if (mode != mode_flag) {
					FATAL ("Loading simulation in different mode.");
					throw 0;
				}
			}
			
			template class boussinesq_element <double>;
		} /* cosine */
	} /* fourier */
} /* two_d */