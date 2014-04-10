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
			boussinesq_element <datatype>::boussinesq_element (bases::axis i_axis_n, bases::axis i_axis_m, int i_name, io::parameters& i_params, bases::messenger* i_messenger_ptr, int i_element_flags) : 
			element <datatype> (i_axis_n, i_axis_m, i_name, i_params, i_messenger_ptr, i_element_flags) {
		
				assert (n > 0);
				assert (m > 0);
		
				TRACE ("Initializing...");
				x_ptr = ptr (x_position);
				z_ptr = ptr (z_position);
				initialize (temp, "T");
				x_vel_ptr = initialize (x_velocity, "u");
				z_vel_ptr = initialize (z_velocity, "w");
				initialize (pressure, "P");
				
				advection_coeff = std::max (i_params.get <datatype> ("temperature.advection"), i_params.get <datatype> ("velocity.advection"));
				cfl = i_params.get <datatype> ("time.cfl");

				DEBUG ("Adding velocity solve.");

				// Solve velocity
				element <datatype>::add_solver (x_velocity, new divergence_solver <datatype> (*grids [0], *grids [1], ptr (x_velocity), ptr (z_velocity), &element_flags [state], &element_flags [x_velocity]));
				
				DEBUG ("Adding z velocity solve.");
				
				element <datatype>::add_solver (z_velocity, new solver <datatype> (*grids [0], *grids [1], messenger_ptr, timestep, alpha_0, alpha_n, ptr (z_velocity), &element_flags [state], &element_flags [z_velocity]));
				
				DEBUG ("Adding velocity plans.");
				
				
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
				
				DEBUG ("Added velocity solve");
				
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
			datatype boussinesq_element <datatype>::calculate_timestep (int i, int j, io::virtual_dump *dump) {
				if (dump) {
					if (j == 0 || j == dump->dims ["z"] [1] - 1) {
						return 1.0 / 0.0;
					}
					if (i == 0 || i == dump->dims ["z"] [0] - 1) {
						return std::abs ((dump->index <datatype> ("z", i, j + 1) - dump->index <datatype> ("z", i, j - 1)) / dump->index <datatype> ("w", i, j) / advection_coeff) * cfl;
					} else {
						return std::min (std::abs ((dump->index <datatype> ("x", i + 1, j) - dump->index <datatype> ("x", i - 1, j)) / dump->index <datatype> ("u", i, j) / advection_coeff), std::abs ((dump->index <datatype> ("z", i, j + 1) - dump->index <datatype> ("z", i, j - 1)) / dump->index <datatype> ("w", i, j) / advection_coeff)) * cfl;
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
		} /* cosine */
	} /* fourier */
} /* two_d */