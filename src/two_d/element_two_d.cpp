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
				
				initialize (temp, &init [0]);
				for (int i = 0; i < n; ++i) {
					init [i] = scale * ((*this) (x_position, i, 0)) + 1.0;
				}
				initialize (x_velocity, &init [0], uniform_n);
				initialize (z_velocity, &init [0], uniform_n & uniform_m);
				initialize (temp_explicit_rhs, NULL, no_transform);
				initialize (temp_implicit_rhs, NULL, no_transform);
				initialize (temp_real_rhs, NULL, only_forward_horizontal);

			
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
				normal_stream->template append <datatype> ("T", ptr (temp));
			
				filestream.str ("");
				filestream << "../output/" + params.output + "_t_" << std::setfill ('0') << std::setw (2) << name << "_%04i.cdf";
				transform_stream.reset (new io::incremental (new io::two_d::netcdf (n, m), filestream.str (), params.output_every));
				transform_stream->template append <int> ("i", &cell_n [0]);
				transform_stream->template append <int> ("j", &cell_m [0]);
				transform_stream->template append <datatype> ("x", ptr (x_position));
				transform_stream->template append <datatype> ("z", ptr (z_position));
				transform_stream->template append <datatype> ("u", ptr (x_velocity));
				transform_stream->template append <datatype> ("w", ptr (z_velocity));
				transform_stream->template append <datatype> ("T", ptr (temp));
							
				// Set up solver
				element <datatype>::add_solver (temp, new solver <datatype> (*grids [0], *grids [1], messenger_ptr, timestep, alpha_0, alpha_n, ptr (temp), ptr (temp_explicit_rhs), ptr (temp_real_rhs), ptr (temp_implicit_rhs)));
		
				// Set up plans in order
				element <datatype>::add_pre_plan (new vertical_diffusion <datatype> (*grids [0], *grids [1], params.diffusion_coeff, params.implicit_alpha, matrix_ptr (temp, 0), matrix_ptr (temp, 1), ptr (temp), ptr (temp_implicit_rhs)));
				element <datatype>::add_mid_plan (new horizontal_diffusion <datatype> (*grids [0], *grids [1], params.diffusion_coeff, params.implicit_alpha, matrix_ptr (temp, 0), matrix_ptr (temp, 1), ptr (temp), ptr (temp_implicit_rhs)));
				if (params.advection_coeff != 0.0) {
					element <datatype>::add_post_plan (new advection <datatype> (*grids [0], *grids [1], params.advection_coeff, ptr (x_velocity), ptr (z_velocity), ptr (temp), ptr (temp_real_rhs)));
				}
		
				TRACE ("Initialized.");
			}
			
			template <class datatype>
			datatype advection_diffusion_element <datatype>::calculate_timestep () {
				datatype t_timestep = params.max_timestep;
				datatype *x_ptr = &((*this) (x_position)), *z_ptr = &((*this) (z_position)), *temp_ptr = &((*this) (temp)), *x_vel_ptr = &((*this) (x_vel)), *z_vel_ptr = &((*this) (z_vel));
				if (params.advection_coeff != 0.0) {
					for (int i = 1; i < n - 1; ++i) {
						for (int j = 1; j < m - 1; ++j) {
							t_timestep = std::min (t_timestep, (datatype) (std::abs ((x_ptr [(i + 1) * m + j] - x_ptr [(i - 1) * m + j]) / x_vel_ptr [i * m + j] / params.advection_coeff) * params.courant_factor));
							t_timestep = std::min (t_timestep, (datatype) (std::abs ((z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / z_vel_ptr [i * m + j] / params.advection_coeff) * params.courant_factor));
						}
					}
				}
				if (params.diffusion_coeff != 0.0) {
					for (int i = 1; i < n - 1; ++i) {
						for (int j = 1; j < m - 1; ++j) {
							t_timestep = std::min (t_timestep, (datatype) (std::abs ((temp_ptr [i * m + j]) / (temp_ptr [(i + 1) * m + j] - 2.0 * temp_ptr [i * m + j] + temp_ptr [(i - 1) * m + j]) * (x_ptr [(i + 1) * m + j] - x_ptr [(i - 1) * m + j]) * (x_ptr [(i + 1) * m + j] - x_ptr [(i - 1) * m + j]) / params.diffusion_coeff)) * params.implicit_allowance);
							t_timestep = std::min (t_timestep, (datatype) (std::abs ((temp_ptr [i * m + j]) / (temp_ptr [i * m + j + 1] - 2.0 * temp_ptr [i * m + j] + temp_ptr [i * m + j - 1]) * (z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) * (z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / params.diffusion_coeff)) * params.implicit_allowance);
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
			
			template <class datatype>
			convection_element <datatype>::convection_element (bases::axis *i_axis_n, bases::axis *i_axis_m, int i_name, io::parameters <datatype>& params, bases::messenger* i_messenger_ptr, int i_flags) : 
			element <datatype> (i_axis_n, i_axis_m, i_name, params, i_messenger_ptr, i_flags) {
		
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
				initialize (temp, &init [0]);
				for (int i = 0; i < n; ++i) {
					init [i] = scale * ((*this) (x_position, i, 0)) + 1.0;
				}
				initialize (stream, &init [0], uniform_n);
				initialize (stream_explicit_rhs, NULL, no_transform);
				initialize (stream_implicit_rhs, NULL, no_transform);
				initialize (stream_real_rhs, NULL, only_forward_horizontal);
				initialize (temp_explicit_rhs, NULL, no_transform);
				initialize (temp_implicit_rhs, NULL, no_transform);
				initialize (temp_real_rhs, NULL, only_forward_horizontal);
			
				// Set up output
				std::ostringstream filestream;
				filestream << "../output/" + params.output + "_" << std::setfill ('0') << std::setw (2) << name << "_%04i.cdf";
				normal_stream.reset (new io::incremental (new io::two_d::netcdf (n, m), filestream.str (), params.output_every));
				normal_stream->template append <int> ("i", &cell_n [0]);
				normal_stream->template append <int> ("j", &cell_m [0]);
				normal_stream->template append <datatype> ("x", ptr (x_position));
				normal_stream->template append <datatype> ("z", ptr (z_position));
				normal_stream->template append <datatype> ("phi", ptr (stream));
				normal_stream->template append <datatype> ("T", ptr (temp));
			
				// Set up solver
				element <datatype>::add_solver (temp, new solver <datatype> (*grids [0], *grids [1], messenger_ptr, timestep, alpha_0, alpha_n, ptr (temp), ptr (temp_explicit_rhs), ptr (temp_real_rhs), ptr (temp_implicit_rhs)));
		
				// Set up plans in order
				element <datatype>::add_pre_plan (new vertical_diffusion <datatype> (*grids [0], *grids [1], params.diffusion_coeff, params.implicit_alpha, matrix_ptr (temp, 0), matrix_ptr (temp, 1), ptr (temp), ptr (temp_implicit_rhs)));
				element <datatype>::add_mid_plan (new horizontal_diffusion <datatype> (*grids [0], *grids [1], params.diffusion_coeff, params.implicit_alpha, matrix_ptr (temp, 0), matrix_ptr (temp, 1), ptr (temp), ptr (temp_implicit_rhs)));
				if (params.advection_coeff != 0.0) {
					element <datatype>::add_post_plan (new stream_advection <datatype> (*grids [0], *grids [1], params.advection_coeff, ptr (stream), ptr (temp), ptr (temp_real_rhs)));
				}
		
				TRACE ("Initialized.");
			}
			
			template <class datatype>
			datatype convection_element <datatype>::calculate_timestep () {
				datatype t_timestep = params.max_timestep;
				datatype *x_ptr = &((*this) (x_position)), *z_ptr = &((*this) (z_position)), *temp_ptr = &((*this) (temp)), *stream_ptr = &((*this) (stream));
				if (params.advection_coeff != 0.0) {
					for (int i = 1; i < n - 1; ++i) {
						for (int j = 1; j < m - 1; ++j) {
							t_timestep = std::min (t_timestep, (datatype) (std::abs ((x_ptr [(i + 1) * m + j] - x_ptr [(i - 1) * m + j]) / (stream_ptr [i * m + j + 1] - stream_ptr [i * m + j - 1]) * (z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / params.advection_coeff) * params.courant_factor));
							t_timestep = std::min (t_timestep, (datatype) (std::abs ((z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / (stream_ptr [(i + 1) * m + j] - stream_ptr [(i - 1) * m + j]) * (x_ptr [(i + 1) * m + j] - x_ptr [(i - 1) * m + j]) / params.advection_coeff) * params.courant_factor));
						}
					}
				}
				if (params.diffusion_coeff != 0.0) {
					for (int i = 1; i < n - 1; ++i) {
						for (int j = 1; j < m - 1; ++j) {
							t_timestep = std::min (t_timestep, (datatype) (std::abs ((temp_ptr [i * m + j]) / (temp_ptr [(i + 1) * m + j] - 2.0 * temp_ptr [i * m + j] + temp_ptr [(i - 1) * m + j]) * (x_ptr [(i + 1) * m + j] - x_ptr [(i - 1) * m + j]) * (x_ptr [(i + 1) * m + j] - x_ptr [(i - 1) * m + j]) / params.diffusion_coeff)) * params.implicit_allowance);
							t_timestep = std::min (t_timestep, (datatype) (std::abs ((temp_ptr [i * m + j]) / (temp_ptr [i * m + j + 1] - 2.0 * temp_ptr [i * m + j] + temp_ptr [i * m + j - 1]) * (z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) * (z_ptr [i * m + j + 1] - z_ptr [i * m + j - 1]) / params.diffusion_coeff)) * params.implicit_allowance);
						}
					}
				}
				if (t_timestep < timestep || t_timestep > 2.0 * timestep) {
					return t_timestep;
				} else {
					return timestep;
				}
			}
			
			template class convection_element <double>;
		} /* chebyshev */
	} /* fourier */
} /* two_d */