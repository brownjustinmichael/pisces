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
			
				// Set up output
				std::ostringstream convert;
				convert << name;
				normal_stream.reset (new io::incremental_output <datatype>  ("../output/normal_" + convert.str () + "_", ".dat", 4, new io::header, n, inputParams["output_every"].asInt));
				normal_stream->append (cell_n [0]);
				normal_stream->append (cell_m [0]);
				normal_stream->append ((*this) [x_position]);
				normal_stream->append ((*this) [z_position]);
				normal_stream->append ((*this) [velocity]);
				normal_stream->append ((*this) [vel_explicit_rhs]);
			
				// Set up plans in order
				element <datatype>::add_pre_plan (new diffusion <datatype> (*grids [0], *grids [1], diffusion_coeff, alpha, pointer (velocity), pointer (vel_implicit_rhs)));

				element <datatype>::add_transform (new transform <datatype> (*grids [0], *grids [1], pointer (velocity)));
		
				// Set up solver
				element <datatype>::add_solver (new solver <datatype> (*grids [0], *grids [1], messenger_ptr, inputParams["n_iterations"].asInt, timestep, pointer (velocity), pointer (vel_explicit_rhs), pointer (vel_implicit_rhs)));
		
				normal_stream->to_file ();
		
				TRACE ("Initialized.");
			}
			
			template <class datatype>
			datatype advection_diffusion_element <datatype>::calculate_timestep () {
				return 0.0;
			}
			
			template class element <float>;
			template class element <double>;
			
			template class advection_diffusion_element <float>;
			template class advection_diffusion_element <double>;
		} /* chebyshev */
	} /* fourier */
} /* two_d */