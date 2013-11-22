/*!***********************************************************************
 * \file element_one_d.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cmath>
#include <cassert>
#include <string>
#include <algorithm>
#include <memory>
#include "../config.hpp"
#include "element_one_d.hpp"
#include "diffusion_one_d.hpp"
#include "advection_one_d.hpp"
#include "solver_one_d.hpp"
#include "transform_one_d.hpp"
	
namespace one_d
{
	namespace chebyshev
	{
		template <class datatype>
		advection_diffusion_element <datatype>::advection_diffusion_element (bases::axis *i_axis_n, int i_name, io::parameters <datatype>& params, bases::messenger* i_messenger_ptr, int i_flags) : 
		element <datatype> (i_axis_n, i_name, params, i_messenger_ptr, i_flags) {
			datatype diffusion_coeff = params.diffusion_coeff;
			datatype advection_coeff = params.advection_coeff; 
			datatype alpha = 0.5;
		
			assert (n > 0);
		
			TRACE ("Initializing...");
			
			datatype scale = params.scale;
			datatype width = params.width;
			datatype mean = params.mean;
			datatype sigma = params.sigma;
			std::vector <datatype> init (n);
			datatype height, temp;
			height = std::max (scale * std::exp (- (width / 2.0 - mean) * (width / 2.0 - mean) / 2.0 / sigma / sigma), scale * std::exp (- (- width / 2.0 - mean) * (- width / 2.0 - mean) / 2.0 / sigma / sigma));
			for (int i = 0; i < n; ++i) {
				temp = scale * std::exp (- ((*this) (position, i) - mean) * ((*this) (position, i) - mean) / 2.0 / sigma / sigma) - height;
				if (temp > 0.0) {
					init [i] = temp;
				} else {
					init [i] = 0.0;
				}
			}
			initialize (velocity, &init [0]);
			initialize (vel_explicit_rhs);
			initialize (vel_implicit_rhs);
			
			// Set up output
			std::ostringstream convert;
			convert << name;
			normal_stream.reset (new io::incremental (new io::one_d::ascii (n), "../output/normal_%04i.dat", params.output_every));
			normal_stream->template append <int> ("i", &(cell [0]));
			normal_stream->template append <datatype> ("x", pointer (position));
			normal_stream->template append <datatype> ("u", pointer (velocity));
			normal_stream->template append <datatype> ("rhs", pointer (vel_explicit_rhs));
			
			// Set up plans in order
			element <datatype>::add_pre_plan (new diffusion <datatype> (*grids [0], diffusion_coeff, alpha, pointer (velocity), pointer (vel_implicit_rhs)));

			element <datatype>::add_inverse_vertical_transform (new fftw_cosine <datatype> (*grids [0], pointer (velocity)));
			element <datatype>::add_forward_vertical_transform (new fftw_cosine <datatype> (*grids [0], pointer (velocity)));
			if (advection_coeff != 0.0) {
				element <datatype>::add_post_plan (new advec <datatype> (*grids [0], advection_coeff, pointer (velocity), pointer (vel_explicit_rhs)));
			}
		
			// Set up solver
			element <datatype>::add_solver (new solver <datatype> (*grids [0], messenger_ptr, timestep, pointer (velocity), pointer (vel_explicit_rhs), pointer (vel_implicit_rhs)));
		
			normal_stream->to_file ();
		
			TRACE ("Initialized.");
		}
		
		template <class datatype>
		datatype advection_diffusion_element <datatype>::calculate_timestep () {
			datatype t_timestep;
			t_timestep = params.max_timestep;
			for (int i = 1; i < n - 1; ++i) {
				t_timestep = std::min (t_timestep, (datatype) (std::abs (((*this) (position, i - 1) - (*this) (position, i + 1)) / (*this) (velocity, i)) / params.advection_coeff));
			}
			t_timestep *= params.courant_factor;
			if (t_timestep < timestep || t_timestep > 2.0 * timestep) {
				return t_timestep;
			} else {
				return timestep;
			}
		}
		
		template class advection_diffusion_element <double>;
	} /* chebyshev */
} /* one_d */