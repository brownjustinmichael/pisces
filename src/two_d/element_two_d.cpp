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
				initialize (velocity, &init [0]);
				initialize (vel_explicit_rhs, NULL, no_transform);
				initialize (vel_implicit_rhs, NULL, no_transform);
			
				// Set up output
				std::ostringstream filestream;
				filestream << "../output/" + params.output + "_" << std::setfill ('0') << std::setw (2) << name << "_%04i.cdf";
				normal_stream.reset (new io::incremental (new io::two_d::netcdf (n, m), filestream.str (), params.output_every));
				normal_stream->template append <int> ("i", &cell_n [0]);
				normal_stream->template append <int> ("j", &cell_m [0]);
				normal_stream->template append <datatype> ("x", pointer (x_position));
				normal_stream->template append <datatype> ("z", pointer (z_position));
				normal_stream->template append <datatype> ("w", pointer (velocity));
				normal_stream->template append <datatype> ("rhs", pointer (vel_implicit_rhs));
			
				// Set up plans in order
				element <datatype>::add_pre_plan (new vertical_diffusion <datatype> (*grids [0], *grids [1], z_diffusion_coeff, alpha, pointer (velocity), pointer (vel_implicit_rhs)));
				element <datatype>::add_mid_plan (new horizontal_diffusion <datatype> (*grids [0], *grids [1], x_diffusion_coeff, alpha, pointer (velocity), pointer (vel_implicit_rhs)));
		
				// Set up solver
				element <datatype>::add_solver (new solver <datatype> (*grids [0], *grids [1], messenger_ptr, timestep, alpha_0, alpha_n, pointer (velocity), pointer (vel_explicit_rhs), pointer (vel_implicit_rhs)));
		
				normal_stream->to_file ();
				
				flags |= x_solve;
		
				TRACE ("Initialized.");
			}
			
			template <class datatype>
			datatype advection_diffusion_element <datatype>::calculate_timestep () {
				return params.max_timestep;
			}
			
			template class element <double>;
			
			template class advection_diffusion_element <double>;
		} /* chebyshev */
	} /* fourier */
} /* two_d */