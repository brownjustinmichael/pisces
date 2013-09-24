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
#include "advection_one_d.h"
#include "solver_one_d.hpp"
#include "fftw_one_d.hpp"
	
namespace one_d
{
	namespace chebyshev
	{
		template <class datatype>
		advection_diffusion_element <datatype>::advection_diffusion_element (int i_n, int i_excess_0, datatype i_position_0, int i_excess_n, datatype i_position_n, int i_name, io::parameter_map& inputParams, bases::messenger* i_messenger_ptr, int i_flags) : 
		element <datatype> (i_n, i_excess_0, i_position_0, i_excess_n, i_position_n, i_name, inputParams, i_messenger_ptr, i_flags) {
			datatype diffusion_coeff = inputParams["diffusion_coeff"].asDouble;
			datatype advection_coeff = inputParams["advection_coeff"].asDouble; 
			datatype alpha = 0.5;
		
			assert (n > 0);
		
			TRACE ("Initializing...");
		
			matrix.resize (n * n, 0.0);
			
			// Set up output
			std::ostringstream convert;
			convert << name;
			normal_stream.reset (new io::incremental_output <datatype>  ("../output/normal_" + convert.str () + "_", ".dat", 4, new io::header, n, inputParams["output_every"].asInt));
			normal_stream->append (cell [0]);
			normal_stream->append ((*this) [position]);
			normal_stream->append ((*this) [velocity]);
			normal_stream->append ((*this) [vel_explicit_rhs]);
			
			// Set up plans in order
			element <datatype>::add_pre_plan (new explicit_diffusion <datatype> (diffusion_coeff * (1.0 - alpha), n, &*grid, pointer (velocity), pointer (vel_implicit_rhs)));

			element <datatype>::add_transform (new fftw_cosine <datatype> (n, pointer (velocity)));
			if (advection_coeff != 0.0) {
				element <datatype>::add_post_plan (new advec <datatype> (n, advection_coeff, pointer (velocity), pointer (vel_explicit_rhs), grid));
			}
			element <datatype>::add_implicit_plan (new implicit_diffusion <datatype> (- diffusion_coeff * alpha, n, &*grid, &matrix [0]));
		
			/*
				TODO Can we include the explicit part of diffusion in the implicit diffusion plan?
			*/
		
			// Set up solver
			element <datatype>::add_solver (new solver <datatype> (messenger_ptr, n, i_excess_0, i_excess_n, inputParams["n_iterations"].asInt, timestep, boundary_weights [edge_0], boundary_weights [edge_n], pointer (position), grid->get_data (0), &matrix [0], pointer (velocity), pointer (vel_explicit_rhs), pointer (vel_implicit_rhs)));
		
			normal_stream->to_file ();
		
			TRACE ("Initialized.");
		}
		
		template <class datatype>
		datatype advection_diffusion_element <datatype>::calculate_timestep () {
			datatype t_timestep;
			t_timestep = inputParams["time_step_size"].asDouble;
			for (int i = 1; i < n - 1; ++i) {
				t_timestep = std::min (t_timestep, (datatype) (std::abs (((*this) (position, i - 1) - (*this) (position, i + 1)) / (*this) (velocity, i)) / inputParams["advection_coeff"].asDouble));
			}
			t_timestep *= inputParams["courant_factor"].asDouble;
			if (t_timestep < timestep || t_timestep > 2.0 * timestep) {
				return t_timestep;
			} else {
				return timestep;
			}
		}
		
		template class advection_diffusion_element <double>;
		template class advection_diffusion_element <float>;
	} /* chebyshev */
} /* one_d */