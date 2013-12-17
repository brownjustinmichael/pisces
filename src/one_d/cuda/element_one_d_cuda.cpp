/*!**********************************************************************
 * \file element_one_d_cuda.cpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-21.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "element_one_d_cuda.hpp"
#include "transform_one_d_cuda.hpp"
#include "solver_one_d_cuda.hpp"
#include "diffusion_one_d_cuda.hpp"
#include "../diffusion_one_d.hpp"

namespace cuda
{
	namespace one_d
	{
		namespace chebyshev
		{
			template <class datatype>
			fft_element <datatype>::fft_element (int i_n, int i_excess_0, datatype i_position_0, int i_excess_n, datatype i_position_n, int i_name, io::parameters <datatype>& params, bases::messenger* i_messenger_ptr, int i_flags) : 
			::one_d::chebyshev::element <datatype> (i_n, i_excess_0, i_position_0, i_excess_n, i_position_n, i_name, params, i_messenger_ptr, i_flags),
			excess_0 (i_excess_0),
			excess_n (i_excess_n) {

				assert (n > 0);
				
				TRACE ("Initializing...");
				
				matrix.resize (n * n, 0.0);
			
				// data_dev is twice necessary size for FFT
				data_dev.resize (2 * n);
				rhs_dev.resize (n);
				
				data_dev.copy_to_device (n, ptr (velocity));
		
				TRACE ("Initialized.");
			}
	
			template <class datatype>
			void fft_element <datatype>::setup () {
				datatype diffusion_coeff = (datatype) params["diffusion_coeff"].asDouble;
				datatype alpha = 0.5;
				
				// Set up output
				std::ostringstream convert;
				convert << name;
				normal_stream.reset (new io::incremental_output <datatype>  ("../output/normal_" + convert.str () + "_", ".dat", 4, new io::header, n, params["output_every"].asInt));
				normal_stream->template append <int> ("i", &cell [0]);
				normal_stream->template append <datatype> ("x", ptr (position));
				normal_stream->template append <datatype> ("u", ptr (velocity));
				
				transform_stream.reset (new io::incremental_output <datatype>  ("../output/transform_" + convert.str () + "_", ".dat", 4, new io::header, n, params["output_every"].asInt));
				normal_stream->template append <int> ("i", &cell [0]);
				normal_stream->template append <datatype> ("x", ptr (position));
				normal_stream->template append <datatype> ("u", ptr (velocity));
		
				transform_stream->to_file ();
				
				this->add_implicit_plan (new ::one_d::implicit_diffusion <datatype> (- diffusion_coeff * alpha, n, &*grid, &matrix [0]));
				
				this->add_pre_plan (new explicit_diffusion <datatype> (diffusion_coeff * (1.0 - alpha), n, &*grid, data_dev.ptr (), rhs_dev.ptr ()));
		
				this->add_transform (new fftw_cosine <datatype> (n, data_dev.ptr ()));
				this->add_post_plan (new transfer <datatype> (n, data_dev.ptr (), &((*this) [velocity])));
				
				// Set up solver
				this->add_solver (new solver <datatype> (messenger_ptr, n, excess_0, excess_n, timestep, boundary_weights [::one_d::edge_0], boundary_weights [::one_d::edge_n], ptr (position), grid->get_data (0), &matrix [0], data_dev.ptr (), rhs_dev.ptr ()));
		
			}
	
			template <class datatype>
			datatype fft_element <datatype>::calculate_timestep () {
				datatype t_timestep;
				t_timestep = params["time_step_size"].asDouble;
				for (int i = 1; i < n - 1; ++i) {
					t_timestep = std::min (t_timestep, (datatype) (std::abs (((*this) (position, i - 1) - (*this) (position, i + 1)) / (*this) (velocity, i))));
				}
				t_timestep *= params["courant_factor"].asDouble;
				if (t_timestep < timestep || t_timestep > 2.0 * timestep) {
					return t_timestep;
				} else {
					return timestep;
				}
			}

			template class fft_element <double>;
			template class fft_element <float>;
		} /* chebyshev */
	} /* one_d */
} /* cuda */
