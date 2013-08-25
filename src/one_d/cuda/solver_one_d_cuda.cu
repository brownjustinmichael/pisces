/*!**********************************************************************
 * \file solver_one_d_cuda.cu
 * /Users/justinbrown/Dropbox/spectral_element/src
 * 
 * Created by Justin Brown on 2013-08-21.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "solver_one_d_cuda.hpp"
#include "fftw_one_d_cuda.hpp"
#include "../../utils/cuda/solver_utils_cuda.hpp"

namespace one_d
{
	namespace cuda
	{
		template <class datatype>
		solver <datatype>::solver (bases::element <datatype>* i_element_ptr, int i_n, int i_excess_0, int i_excess_n, datatype& i_timestep, datatype& i_alpha_0, datatype& i_alpha_n, datatype *i_default_matrix, datatype *i_matrix, int i_name_in, int i_name_rhs, int i_name_out = null, int i_flags = 0x00) :
		one_d::solver <datatype> (i_element_ptr, i_n, i_excess_0, i_excess_n, i_timestep, i_alpha_0, i_alpha_n, i_default_matrix, i_matrix, i_name_in, i_name_rhs, i_name_out, i_flags) {
			HANDLE_ERROR (cudaMalloc ((void**) &factorized_matrix_dev, n * n * sizeof (datatype)));
			HANDLE_ERROR (cudaMalloc ((void**) &ipiv_dev, n * sizeof (int)));
		}
		
		template <class datatype>
		solver <datatype>::~solver () {
			HANDLE_ERROR (cudaFree (factorized_matrix_dev));
			HANDLE_ERROR (cudaFree (ipiv_dev));
		}
		
		template <class datatype>
		void solver <datatype>::_factorize () {
			one_d::solver <datatype>::_factorize ();
			
			HANDLE_ERROR (cudaMemcpy (factorized_matrix_dev, &factorized_matrix [0], n * n * sizeof (datatype), cudaMemcpyHostToDevice));
			HANDLE_ERROR (cudaMemcpy (ipiv_dev, &ipiv [0], n * sizeof (int), cudaMemcpyHostToDevice));
		}
		
		template <class datatype>
		void solver <datatype>::execute () {
			utils::cuda::copy (n, data_in, data_out);
		
			data_out [0] += alpha_0 * timestep * rhs [0];
			data_out [n - 1] += alpha_n * timestep * rhs [n - 1];
			utils::add_scaled (n - 2, timestep, rhs + 1, data_out + 1);
			
			utils::cuda::matrix_solve (n, factorized_matrix_dev, ipiv_dev, data_out);
		}
		
		template class solver <double>;
	} /* cuda */
} /* one_d */