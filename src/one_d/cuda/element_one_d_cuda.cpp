/*!**********************************************************************
 * \file element_one_d_cuda.cpp
 * /Users/justinbrown/Dropbox/spectral_element/src
 * 
 * Created by Justin Brown on 2013-08-21.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "element_one_d_cuda.hpp"
#include "fftw_one_d_cuda.hpp"
#include "solver_one_d_cuda.hpp"

namespace one_d
{
	namespace chebyshev
	{
		template <class datatype>
		class element;
		
		namespace cuda
		{
			template <class datatype>
			fft_element <datatype>::fft_element (int i_n, int i_excess_0, datatype i_position_0, int i_excess_n, datatype i_position_n, int i_name, io::parameter_map& inputParams, bases::messenger <datatype>* i_messenger_ptr, int i_flags) : 
			element <datatype> (i_n, i_excess_0, i_position_0, i_excess_n, i_position_n, i_name, inputParams, i_messenger_ptr, i_flags),
			excess_0 (i_excess_0),
			excess_n (i_excess_n) {

				assert (n > 0);
					
				TRACE ("Initializing...");
					
				matrix.resize (i_n * i_n, 0.0);
				
				// data_dev is twice necessary size for FFT
				data_dev.resize (2 * n);
				rhs_dev.resize (n);
			
				TRACE ("Initialized.");
			}
		
			template <class datatype>
			void fft_element <datatype>::setup () {
				// Set up output
				std::ostringstream convert;
				convert << name;
				transform_stream.reset (new io::incremental_output <datatype>  ("../output/test_angle_" + convert.str () + "_", ".dat", 4, new io::header, n, inputParams["output_every"].asInt));
				transform_stream->append (cell [0]);
				transform_stream->append ((*this) [position]);
				transform_stream->append ((*this) [velocity]);
				transform_stream->append ((*this) [rhs]);
			
				transform_stream->to_file ();
			
				element <datatype>::add_transform (new one_d::cuda::fftw_cosine (this, n, data_dev.pointer ()));
				element <datatype>::add_post_plan (new one_d::cuda::transfer (this, n, data_dev.pointer (), &((*this) [velocity])));
					
				// Set up solver
				element <datatype>::add_solver (new solver <datatype> (this, n, excess_0, excess_n, timestep, boundary_weights [edge_0], boundary_weights [edge_n], grid->get_data (0), &matrix [0], velocity, rhs));
			
			}
		
			template <class datatype>
			datatype fft_element <datatype>::calculate_timestep () {
				return 0.0;
			}

			template class fft_element <double>;
		} /* cuda */
	} /* chebyshev */
} /* one_d */