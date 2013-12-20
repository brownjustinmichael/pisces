/*!***********************************************************************
 * \file element_one_d.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "element_one_d.hpp"

#include <cmath>
#include <cassert>
#include <string>
#include <algorithm>
#include <memory>
#include <iomanip>
#include "../config.hpp"
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
			datatype alpha = params.implicit_alpha;
		
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
			initialize (vel_explicit_rhs, NULL, no_transform);
			initialize (vel_implicit_rhs, NULL, no_transform);
			
			// Set up output
			std::ostringstream filestream;
			filestream << "../output/output_" << std::setfill ('0') << std::setw (2) << name << "_%04i.dat";
			normal_stream.reset (new io::incremental (new io::one_d::ascii (n), filestream.str (), params.output_every));
			normal_stream->template append <int> ("i", &(cell [0]));
			normal_stream->template append <datatype> ("x", ptr (position));
			normal_stream->template append <datatype> ("u", ptr (velocity));
			normal_stream->template append <datatype> ("rhs", ptr (vel_explicit_rhs));
			
			// Set up solver
			element <datatype>::add_solver (velocity, new solver <datatype> (*grids [0], messenger_ptr, timestep, alpha_0, alpha_n, ptr (velocity), ptr (vel_explicit_rhs), ptr (vel_implicit_rhs)));
			
			// Set up plans in order
			element <datatype>::add_pre_plan (new diffusion <datatype> (*grids [0], diffusion_coeff, alpha, matrix_ptr (velocity), ptr (velocity), ptr (vel_implicit_rhs)));
			if (advection_coeff != 0.0) {
				element <datatype>::add_post_plan (new advec <datatype> (*grids [0], advection_coeff, ptr (velocity), ptr (vel_explicit_rhs)));
			}
			
			normal_stream->to_file ();
		
			TRACE ("Initialized.");
		}
		
		template <class datatype>
		datatype advection_diffusion_element <datatype>::calculate_timestep () {
			datatype t_timestep;
			t_timestep = params.max_timestep;
			for (int i = 1; i < n - 1; ++i) {
				t_timestep = std::min (t_timestep, (datatype) (std::abs (((*this) (position, i - 1) - (*this) (position, i + 1)) / (*this) (velocity, i) / params.advection_coeff) * params.courant_factor));
			}
			if (t_timestep < timestep || t_timestep > 2.0 * timestep) {
				return t_timestep;
			} else {
				return timestep;
			}
		}
		
		template class advection_diffusion_element <double>;
		
		template <class datatype>
		nonlinear_diffusion_element <datatype>::nonlinear_diffusion_element (bases::axis *i_axis_n, int i_name, io::parameters <datatype>& params, bases::messenger* i_messenger_ptr, int i_flags) : 
		element <datatype> (i_axis_n, i_name, params, i_messenger_ptr, i_flags) {
			datatype diffusion_coeff = params.diffusion_coeff;
			datatype advection_coeff = params.advection_coeff; 
			datatype alpha = params.implicit_alpha;
		
			assert (n > 0);
		
			TRACE ("Initializing...");
			
			datatype scale = params.scale;
			datatype width = params.width;
			datatype mean = params.mean - 0.5;
			datatype sigma = params.sigma;
			std::vector <datatype> init (n);
			datatype height, temp;
			for (int i = 0; i < n; ++i) {
				// if ((*this) (position, i) < 0.0) {
				// 	init [i] = 0.0;
				// } else {
				// 	init [i] = 1.0;
				// }
				init [i] = ((*this) (position, i) - (*this) (position)) * scale;
			}
			height = std::max (scale * std::exp (- (width / 2.0 - mean) * (width / 2.0 - mean) / 2.0 / sigma / sigma), scale * std::exp (- (- width / 2.0 - mean) * (- width / 2.0 - mean) / 2.0 / sigma / sigma));
			for (int i = 0; i < n; ++i) {
				temp = scale * std::exp (- ((*this) (position, i) - mean) * ((*this) (position, i) - mean) / 2.0 / sigma / sigma) - height;
				if (temp > 0.0) {
					init [i] += temp;
				}
			}
			mean = params.mean + 0.5;
			height = std::max (scale * std::exp (- (width / 2.0 - mean) * (width / 2.0 - mean) / 2.0 / sigma / sigma), scale * std::exp (- (- width / 2.0 - mean) * (- width / 2.0 - mean) / 2.0 / sigma / sigma));
			for (int i = 0; i < n; ++i) {
				temp = scale * std::exp (- ((*this) (position, i) - mean) * ((*this) (position, i) - mean) / 2.0 / sigma / sigma) - height;
				if (temp > 0.0) {
					init [i] += temp;
				}
			}
			initialize (velocity, &init [0]);
			initialize (vel_explicit_rhs, NULL, no_transform);
			initialize (vel_implicit_rhs, NULL, no_transform);
			
			// Set up output
			std::ostringstream filestream;
			filestream << "../output/" + params.output + "_" << std::setfill ('0') << std::setw (2) << name << "_%04i.dat";
			normal_stream.reset (new io::incremental (new io::one_d::ascii (n), filestream.str (), params.output_every));
			normal_stream->template append <int> ("i", &(cell [0]));
			normal_stream->template append <datatype> ("x", ptr (position));
			normal_stream->template append <datatype> ("u", ptr (velocity));
			normal_stream->template append <datatype> ("rhs", ptr (vel_explicit_rhs));
			

			// Set up solver
			element <datatype>::add_solver (velocity, new solver <datatype> (*grids [0], messenger_ptr, timestep, alpha_0, alpha_n, ptr (velocity), ptr (vel_explicit_rhs), ptr (vel_implicit_rhs)));
			
			// Set up plans in order
			element <datatype>::add_pre_plan (new diffusion <datatype> (*grids [0], diffusion_coeff, alpha, matrix_ptr (velocity), ptr (velocity), ptr (vel_implicit_rhs)));
			if (params.nonlinear_diffusion_coeff != 0.0) {
				element <datatype>::add_post_plan (new nonlinear_diffusion <datatype> (*grids [0], params.nonlinear_diffusion_coeff, ptr (velocity), ptr (vel_explicit_rhs)));
			}
			if (advection_coeff != 0.0) {
				element <datatype>::add_post_plan (new advec <datatype> (*grids [0], advection_coeff, ptr (velocity), ptr (vel_explicit_rhs)));
			}
		
			normal_stream->to_file ();
		
			TRACE ("Initialized.");
		}
		
		template <class datatype>
		datatype nonlinear_diffusion_element <datatype>::calculate_timestep () {
			datatype t_timestep;
			t_timestep = params.max_timestep;
			// for (int i = 1; i < n - 1; ++i) {
			// 	t_timestep = std::min (t_timestep, (datatype) (std::abs (((*this) (position, i - 1) - (*this) (position, i + 1)) / (*this) (velocity, i) / params.advection_coeff) * params.courant_factor));
			// 	t_timestep = std::min (t_timestep, (datatype) std::abs (((*this) (position, i + 1) - (*this) (position, i - 1)) * ((*this) (position, i + 1) - (*this) (position, i - 1)) / 2.0 / params.nonlinear_diffusion_coeff / (*this) (velocity, i) * params.courant_factor));
			// }
			if (t_timestep < timestep || t_timestep > 2.0 * timestep) {
				return t_timestep;
			} else {
				return timestep;
			}
		}
		
		template class nonlinear_diffusion_element <double>;
	} /* chebyshev */
} /* one_d */