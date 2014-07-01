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
			boussinesq_element <datatype>::boussinesq_element (bases::axis i_axis_n, bases::axis i_axis_m, int i_name, io::parameters& i_params, utils::messenger* i_messenger_ptr, int i_element_flags) : 
			element <datatype> (i_axis_n, i_axis_m, i_name, i_params, i_messenger_ptr, i_element_flags) {
		
				assert (n > 0);
				assert (m > 0);
		
				TRACE ("Initializing...");
				x_ptr = ptr (x_position);
				z_ptr = ptr (z_position);
				initialize (temp, "T");
				initialize (composition, "S");
				x_vel_ptr = initialize (x_velocity, "u");
				z_vel_ptr = initialize (z_velocity, "w");
				initialize (pressure, "P");
				
				advection_coeff = i_params.get <datatype> ("temperature.advection");
				advection_coeff = std::max (advection_coeff, i_params.get <datatype> ("velocity.advection"));
				advection_coeff = std::max (advection_coeff, i_params.get <datatype> ("composition.advection"));
				cfl = i_params.get <datatype> ("time.cfl");
				
				// Solve velocity
				// element <datatype>::add_solver (x_velocity, new divergence_solver <datatype> (*grids [0], *grids [1], ptr (x_velocity), ptr (z_velocity), &element_flags [state], &element_flags [x_velocity]));
				//
				// element <datatype>::add_solver (z_velocity, new solver <datatype> (*grids [0], *grids [1], messenger_ptr, timestep, alpha_0, alpha_n, ptr (z_velocity), &element_flags [state], &element_flags [z_velocity]));
				//
				// solvers [z_velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new vertical_diffusion <datatype> (*solvers [z_velocity], i_params.get <datatype> ("velocity.diffusion"), i_params.get <datatype> ("time.alpha"))), pre_plan);
				// solvers [z_velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new horizontal_diffusion <datatype> (*solvers [z_velocity], i_params.get <datatype> ("velocity.diffusion"), i_params.get <datatype> ("time.alpha"))), mid_plan);
				// solvers [z_velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new advection <datatype> (*solvers [z_velocity], i_params.get <datatype> ("velocity.advection"), ptr (x_velocity), ptr (z_velocity))), post_plan);
				// solvers [z_velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new source <datatype> (*solvers [z_velocity], i_params.get <datatype> ("velocity.buoyancy.temperature"), ptr (temp))), mid_plan);
				// solvers [z_velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new source <datatype> (*solvers [z_velocity], i_params.get <datatype> ("velocity.buoyancy.composition"), ptr (composition))), mid_plan);
				
				// std::shared_ptr <bases::solver <datatype>> pressure_solver = std::shared_ptr <bases::solver <datatype>> (new laplace_solver <datatype> (*grids [0], *grids [1], messenger_ptr, ptr (pressure), &element_flags [state], &element_flags [pressure]));
				// solvers [z_velocity]->add_plan (pressure_solver, pre_solve_plan);
				// solvers [z_velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new z_derivative_source <datatype> (*pressure_solver, 1.0, ptr (temp))), mid_plan);
				// solvers [z_velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new mixed_derivative_source <datatype> (*pressure_solver, -2.0, ptr (z_vel), ptr (x_vel))), post_plan);
				// solvers [z_velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new square_z_derivative_source <datatype> (*pressure_solver, -1.0, ptr (z_vel))), post_plan);
				// solvers [z_velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new square_x_derivative_source <datatype> (*pressure_solver, -1.0, ptr (x_vel))), post_plan);
				// solvers [z_velocity]->add_plan (std::shared_ptr <bases::plan <datatype> > (new z_derivative_source <datatype> (*solvers [z_velocity], -1.0, ptr (pressure))), pre_solve_plan);
				
				// Solve temperature
				element <datatype>::add_solver (temp, new solver <datatype> (*grids [0], *grids [1], messenger_ptr, timestep, alpha_0, alpha_n, ptr (temp), &element_flags [state], &element_flags [temp]));
						
				solvers [temp]->add_plan (std::shared_ptr <bases::plan <datatype>> (new vertical_diffusion <datatype> (*solvers [temp], i_params.get <datatype> ("temperature.diffusion"), i_params.get <datatype> ("time.alpha"))), pre_plan);
				solvers [temp]->add_plan (std::shared_ptr <bases::plan <datatype>> (new horizontal_diffusion <datatype> (*solvers [temp], i_params.get <datatype> ("temperature.diffusion"), i_params.get <datatype> ("time.alpha"))), mid_plan);
				solvers [temp]->add_plan (std::shared_ptr <bases::plan <datatype>> (new advection <datatype> (*solvers [temp], i_params.get <datatype> ("temperature.advection"), ptr (x_vel), ptr (z_vel))), post_plan);
				// solvers [temp]->add_plan (std::shared_ptr <bases::plan <datatype>> (new source <datatype> (*solvers [temp], -i_params.get <datatype> ("temperature.stratification"), ptr (z_velocity))), mid_plan);
				
				// Solve composition
				// element <datatype>::add_solver (composition, new solver <datatype> (*grids [0], *grids [1], messenger_ptr, timestep, alpha_0, alpha_n, ptr (composition), &element_flags [state], &element_flags [composition]));
				//
				// solvers [composition]->add_plan (std::shared_ptr <bases::plan <datatype>> (new vertical_diffusion <datatype> (*solvers [composition], i_params.get <datatype> ("composition.diffusion"), i_params.get <datatype> ("time.alpha"))), pre_plan);
				// solvers [composition]->add_plan (std::shared_ptr <bases::plan <datatype>> (new horizontal_diffusion <datatype> (*solvers [composition], i_params.get <datatype> ("composition.diffusion"), i_params.get <datatype> ("time.alpha"))), mid_plan);
				// solvers [composition]->add_plan (std::shared_ptr <bases::plan <datatype>> (new advection <datatype> (*solvers [composition], i_params.get <datatype> ("composition.advection"), ptr (x_vel), ptr (z_vel))), post_plan);
				// solvers [composition]->add_plan (std::shared_ptr <bases::plan <datatype>> (new source <datatype> (*solvers [composition], -i_params.get <datatype> ("composition.stratification"), ptr (z_velocity))), mid_plan);
				
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