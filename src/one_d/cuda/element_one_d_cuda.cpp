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

namespace cuda
{
	namespace one_d
	{
		namespace chebyshev
		{
			template <class datatype>
			fft_element <datatype>::fft_element (int i_n, int i_excess_0, datatype i_position_0, int i_excess_n, datatype i_position_n, int i_name, io::parameter_map& inputParams, bases::messenger* i_messenger_ptr, int i_flags) : 
			::one_d::chebyshev::element <datatype> (i_n, i_excess_0, i_position_0, i_excess_n, i_position_n, i_name, inputParams, i_messenger_ptr, i_flags),
			excess_0 (i_excess_0),
			excess_n (i_excess_n) {

				assert (n > 0);
				
				TRACE ("Initializing...");
				
				matrix.resize (n * n, 0.0);
			
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
		
				::one_d::chebyshev::element <datatype>::add_transform (new fftw_cosine <datatype> (this, n, data_dev.pointer ()));
				::one_d::chebyshev::element <datatype>::add_post_plan (new transfer <datatype> (this, n, data_dev.pointer (), &((*this) [velocity])));
				
				// Set up solver
				::one_d::chebyshev::element <datatype>::add_solver (new solver <datatype> (this, n, excess_0, excess_n, timestep, boundary_weights [::one_d::edge_0], boundary_weights [::one_d::edge_n], grid->get_data (0), &matrix [0], velocity, rhs));
		
			}
	
			template <class datatype>
			datatype fft_element <datatype>::calculate_timestep () {
				return 0.0;
			}

			template class fft_element <double>;
			template class fft_element <float>;
		} /* chebyshev */
	} /* one_d */
} /* cuda */
