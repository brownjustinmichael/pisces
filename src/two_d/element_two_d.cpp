/*!**********************************************************************
 * \file element_two_d.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-12.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "element_two_d.hpp"
#include "advection_two_d.hpp"
#include "diffusion_two_d.hpp"
#include "transform_two_d.hpp"
#include "solver_two_d.hpp"
#include <sstream>
#include <iomanip>

namespace two_d
{
	namespace fourier
	{
		namespace chebyshev
		{
			template <class datatype>
			advection_diffusion_element <datatype>::advection_diffusion_element (bases::axis *i_axis_n, bases::axis *i_axis_m, int i_name, io::parameters <datatype>& params, bases::messenger* i_messenger_ptr, int i_flags) : 
			element <datatype> (i_axis_n, i_axis_m, i_name, params, i_messenger_ptr, i_flags) {
				datatype x_diffusion_coeff = params.x_diffusion_coeff;
				datatype z_diffusion_coeff = params.z_diffusion_coeff;
				datatype alpha = params.implicit_alpha;
		
				assert (n > 0);
				assert (m > 0);
		
				TRACE ("Initializing...");
				
				datatype scale = params.scale;
				datatype width = params.width;
				datatype mean = params.mean;
				datatype sigma = params.sigma;
				std::vector <datatype> init (n * m);
				datatype height;
				height = std::max (scale * std::exp (- (width / 2.0 - mean) * (width / 2.0 - mean) / 2.0 / sigma / sigma), scale * std::exp (- (- width / 2.0 - mean) * (- width / 2.0 - mean) / 2.0 / sigma / sigma));
				TRACE ("scale " << scale << " width " << width << " mean " << mean << " sigma " << sigma << " height " << height);
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j < m; ++j) {
						init [i * m + j] = scale * (std::exp (- ((*this) (x_position, i, j) - mean) * ((*this) (x_position, i, j) - mean) / 2.0 / sigma / sigma) - height) * (std::exp (- ((*this) (z_position, i, j) - mean) * ((*this) (z_position, i, j) - mean) / 2.0 / sigma / sigma) - height);
					}
				}
				initialize (x_velocity, &init [0]);
				initialize (z_velocity, &init [0]);
				initialize (x_vel_explicit_rhs, NULL, no_transform);
				initialize (x_vel_implicit_rhs, NULL, no_transform);
				initialize (x_vel_real_rhs, NULL, only_forward_horizontal);
				initialize (z_vel_explicit_rhs, NULL, no_transform);
				initialize (z_vel_implicit_rhs, NULL, no_transform);
				initialize (z_vel_real_rhs, NULL, only_forward_horizontal);
			
				// Set up output
				std::ostringstream filestream;
				filestream << "../output/" + params.output + "_" << std::setfill ('0') << std::setw (2) << name << "_%04i.cdf";
				normal_stream.reset (new io::incremental (new io::two_d::netcdf (n, m), filestream.str (), params.output_every));
				normal_stream->template append <int> ("i", &cell_n [0]);
				normal_stream->template append <int> ("j", &cell_m [0]);
				normal_stream->template append <datatype> ("x", ptr (x_position));
				normal_stream->template append <datatype> ("z", ptr (z_position));
				normal_stream->template append <datatype> ("u", ptr (x_velocity));
				normal_stream->template append <datatype> ("w", ptr (z_velocity));
			
				filestream.str ("");
				filestream << "../output/" + params.output + "_t_" << std::setfill ('0') << std::setw (2) << name << "_%04i.cdf";
				transform_stream.reset (new io::incremental (new io::two_d::netcdf (n, m), filestream.str (), params.output_every));
				transform_stream->template append <int> ("i", &cell_n [0]);
				transform_stream->template append <int> ("j", &cell_m [0]);
				transform_stream->template append <datatype> ("x", ptr (x_position));
				transform_stream->template append <datatype> ("z", ptr (z_position));
				transform_stream->template append <datatype> ("u", ptr (x_velocity));
				transform_stream->template append <datatype> ("w", ptr (z_velocity));
							
				// Set up solver
				element <datatype>::add_solver (x_velocity, new solver <datatype> (*grids [0], *grids [1], messenger_ptr, timestep, alpha_0, alpha_n, ptr (x_velocity), ptr (x_vel_explicit_rhs), ptr (x_vel_real_rhs), ptr (x_vel_implicit_rhs)));
				element <datatype>::add_solver (z_velocity, new solver <datatype> (*grids [0], *grids [1], messenger_ptr, timestep, alpha_0, alpha_n, ptr (z_velocity), ptr (z_vel_explicit_rhs), ptr (z_vel_real_rhs), ptr (z_vel_implicit_rhs)));
		
				// Set up plans in order
				element <datatype>::add_pre_plan (new vertical_diffusion <datatype> (*grids [0], *grids [1], x_diffusion_coeff, alpha, matrix_ptr (x_velocity, 0), matrix_ptr (x_velocity, 1), ptr (x_velocity), ptr (x_vel_implicit_rhs)));
				element <datatype>::add_mid_plan (new horizontal_diffusion <datatype> (*grids [0], *grids [1], x_diffusion_coeff, alpha, matrix_ptr (x_velocity, 0), matrix_ptr (x_velocity, 1), ptr (x_velocity), ptr (x_vel_implicit_rhs)));
				if (params.advection_coeff != 0.0) {
					element <datatype>::add_post_plan (new advection <datatype> (*grids [0], *grids [1], params.advection_coeff, ptr (x_velocity), ptr (x_velocity), ptr (x_velocity), ptr (x_vel_real_rhs)));
				}
				
				element <datatype>::add_pre_plan (new vertical_diffusion <datatype> (*grids [0], *grids [1], z_diffusion_coeff, alpha, matrix_ptr (z_velocity, 0), matrix_ptr (z_velocity, 1), ptr (z_velocity), ptr (z_vel_implicit_rhs)));
				element <datatype>::add_mid_plan (new horizontal_diffusion <datatype> (*grids [0], *grids [1], x_diffusion_coeff, alpha, matrix_ptr (z_velocity, 0), matrix_ptr (z_velocity, 1), ptr (z_velocity), ptr (z_vel_implicit_rhs)));
				if (params.advection_coeff != 0.0) {
					element <datatype>::add_post_plan (new advection <datatype> (*grids [0], *grids [1], params.advection_coeff, ptr (x_velocity), ptr (z_velocity), ptr (z_velocity), ptr (z_vel_real_rhs)));
				}
		
				TRACE ("Initialized.");
			}
			
			template <class datatype>
			datatype advection_diffusion_element <datatype>::calculate_timestep () {
				datatype t_timestep = params.max_timestep;
				if (params.advection_coeff != 0.0) {
					for (int i = 1; i < n - 1; ++i) {
						for (int j = 1; j < m - 1; ++j) {
							t_timestep = std::min (t_timestep, (datatype) (0.5 * std::abs (((*this) (x_position, i + 1, j) - (*this) (x_position, i - 1, j)) / (*this) (x_velocity, i, j) / params.advection_coeff) * params.courant_factor));
							t_timestep = std::min (t_timestep, (datatype) (0.5 * std::abs (((*this) (z_position, i, j + 1) - (*this) (z_position, i, j - 1)) / (*this) (z_velocity, i, j) / params.advection_coeff) * params.courant_factor));
							// t_timestep = std::min (t_timestep, (datatype) std::abs (((*this) (position, i + 1) - (*this) (position, i - 1)) * ((*this) (position, i + 1) - (*this) (position, i - 1)) / 2.0 / params.nonlinear_diffusion_coeff / (*this) (velocity, i) * params.courant_factor));
						}
					}
				}
				if (t_timestep < timestep || t_timestep > 2.0 * timestep) {
					return t_timestep;
				} else {
					return timestep;
				}
			}
			
			template class element <double>;
			
			template class advection_diffusion_element <double>;
		} /* chebyshev */
	} /* fourier */
} /* two_d */