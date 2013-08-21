/*!**********************************************************************
 * \file element_one_d_cuda.cpp
 * /Users/justinbrown/Dropbox/spectral_element/src
 * 
 * Created by Justin Brown on 2013-08-21.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "element_one_d_cuda.hpp"
#include "fftw_one_d_cuda.hpp"
#include "../solver_one_d.hpp"

namespace one_d
{
	namespace chebyshev
	{
		template <class datatype>
		cuda_element <datatype>::cuda_element (int i_n, datatype i_position_0, datatype i_position_n, int i_excess_0, int i_excess_n, int i_name, io::parameter_map& inputParams, bases::messenger <datatype>* i_messenger_ptr, int i_flags) : 
		element <datatype> (i_n, return_position (i_n, 0, i_excess_0, i_position_0, i_excess_n, i_position_n), return_position (i_n, i_n - 1, i_excess_0, i_position_0, i_excess_n, i_position_n), i_name, inputParams, i_messenger_ptr, i_flags),
		excess_0 (i_excess_0),
		excess_n (i_excess_n) {

			assert (n > 0);
					
			TRACE ("Initializing...");
					
			matrix.resize (i_n * i_n, 0.0);
			
			TRACE ("Initialized.");

		}
		
		template <class datatype>
		void cuda_element <datatype>::setup () {
			// Set up output
			std::ostringstream convert;
			convert << name;
			normal_stream.reset (new io::incremental_output <datatype>  ("../output/test_angle_" + convert.str () + "_", ".dat", 4, new io::header, n, inputParams["output_every"].asInt));
			normal_stream->append (cell [0]);
			normal_stream->append ((*this) [position]);
			normal_stream->append ((*this) [velocity]);
			normal_stream->append ((*this) [rhs]);
			
			normal_stream->to_file ();
			
			element <datatype>::add_transform (new cuda::fftw_cosine (this, n, velocity));
					
			// Set up solver
			element <datatype>::add_solver (new solver <datatype> (this, n, excess_0, excess_n, timestep, boundary_weights [edge_0], boundary_weights [edge_n], grid->get_data (0), &matrix [0], velocity, rhs));
			
		}
		
		template <class datatype>
		datatype cuda_element <datatype>::calculate_timestep () {
			return 0.0;
		}

		template class cuda_element <double>;
	} /* chebyshev */
} /* one_d */