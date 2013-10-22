/*!**********************************************************************
 * \file element_two_d.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-12.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "element_two_d.hpp"
#include "diffusion_two_d.hpp"
#include "transform_two_d.hpp"
#include "solver_two_d.hpp"

namespace two_d
{
	namespace fourier
	{
		namespace chebyshev
		{
			template <class datatype>
			advection_diffusion_element <datatype>::advection_diffusion_element (bases::axis *i_axis_n, bases::axis *i_axis_m, int i_name, io::parameter_map& inputParams, bases::messenger* i_messenger_ptr, int i_flags) : 
			element <datatype> (i_axis_n, i_axis_m, i_name, inputParams, i_messenger_ptr, i_flags) {
				datatype diffusion_coeff = inputParams["diffusion_coeff"].asDouble;
				datatype alpha = 0.5;
		
				assert (n > 0);
				assert (m > 0);
		
				TRACE ("Initializing...");
				
				datatype scale = inputParams["init_cond_scale"].asDouble;
				datatype width = inputParams["init_cond_width"].asDouble;
				datatype mean = inputParams["init_cond_mean"].asDouble;
				datatype sigma = inputParams["init_cond_sigma"].asDouble;
				std::vector <datatype> init (n * m);
				datatype height;
				height = std::max (scale * std::exp (- (width / 2.0 - mean) * (width / 2.0 - mean) / 2.0 / sigma / sigma), scale * std::exp (- (- width / 2.0 - mean) * (- width / 2.0 - mean) / 2.0 / sigma / sigma));
				DEBUG ("scale " << scale << " width " << width << " mean " << mean << " sigma " << sigma << " height " << height);
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j < m; ++j) {
						init [i * m + j] = scale * std::exp (- (((*this) (x_position, i, j) - mean) * ((*this) (x_position, i, j) - mean) + ((*this) (z_position, i, j) - mean) * ((*this) (z_position, i, j) - mean)) / 2.0 / sigma / sigma) - height;
					}
				}
				initialize (velocity, &init [0]);
				initialize (vel_explicit_rhs);
				initialize (vel_implicit_rhs);
			
				// Set up output
				std::ostringstream convert;
				convert << name;
				normal_stream.reset (new io::incremental (new io::two_d::netcdf (n, m), "../output/normal_%04i.cdf", inputParams["output_every"].asInt));
				normal_stream->template append <int> ("i", &cell_n [0]);
				normal_stream->template append <int> ("j", &cell_m [0]);
				normal_stream->template append <datatype> ("x", pointer (x_position));
				normal_stream->template append <datatype> ("z", pointer (z_position));
				normal_stream->template append <datatype> ("w", pointer (velocity));
			
				transform_stream.reset (new io::incremental (new io::two_d::netcdf (n, m), "../output/transform_%04i.cdf", inputParams["output_every"].asInt));
				transform_stream->template append <int> ("i", &cell_n [0]);
				transform_stream->template append <int> ("j", &cell_m [0]);
				transform_stream->template append <datatype> ("x", pointer (x_position));
				transform_stream->template append <datatype> ("z", pointer (z_position));
				transform_stream->template append <datatype> ("w", pointer (velocity));
				transform_stream->template append <datatype> ("rhs", pointer (vel_implicit_rhs));
			
				// Set up plans in order
				element <datatype>::add_pre_plan (new diffusion <datatype> (*grids [0], *grids [1], diffusion_coeff, alpha, pointer (velocity), pointer (vel_implicit_rhs)));

				element <datatype>::add_forward_transform (new transform <datatype> (*grids [0], *grids [1], pointer (velocity), NULL, ignore_m));
				element <datatype>::add_inverse_transform (new transform <datatype> (*grids [0], *grids [1], pointer (velocity), NULL, ignore_m | inverse));
		
				// Set up solver
				element <datatype>::add_solver (new solver <datatype> (*grids [0], *grids [1], messenger_ptr, inputParams["n_iterations"].asInt, timestep, pointer (velocity), pointer (vel_explicit_rhs), pointer (vel_implicit_rhs)));
		
				normal_stream->to_file ();
		
				TRACE ("Initialized.");
			}
			
			template <class datatype>
			datatype advection_diffusion_element <datatype>::calculate_timestep () {
				return inputParams ["time_step_size"].asDouble;
			}
			
			template class element <float>;
			template class element <double>;
			
			template class advection_diffusion_element <float>;
			template class advection_diffusion_element <double>;
		} /* chebyshev */
	} /* fourier */
} /* two_d */