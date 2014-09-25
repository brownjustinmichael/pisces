/*!**********************************************************************
 * \file solver_one_d_cuda.cu
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-21.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "solver_one_d_cuda.hpp"
#include "transform_one_d_cuda.hpp"
#include "../../utils/cuda/linalg_cuda.hpp"
#include "../../utils/cuda/utils_cuda.cuh"
#include "../../utils/cuda/utils_cublas.hpp"

namespace cuda
{
	namespace one_d
	{
		template <class datatype>
		solver <datatype>::solver (utils::messenger* i_messenger_ptr, int i_n, int i_excess_0, int i_excess_n, datatype& i_timestep, datatype& i_alpha_0, datatype& i_alpha_n, datatype* i_positions, datatype *i_default_matrix, datatype *i_matrix, datatype* i_data_in, datatype* i_rhs, datatype* i_data_out, int i_flags) :
		::one_d::solver <datatype> (i_messenger_ptr, i_n, i_excess_0, i_excess_n, i_timestep, i_alpha_0, i_alpha_n, i_positions, i_default_matrix, i_matrix, i_data_in, i_rhs, i_data_out, i_flags) {
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
			TRACE ("Factorizing...");
			
			::one_d::solver <datatype>::_factorize ();
			
			HANDLE_ERROR (cudaMemcpy (factorized_matrix_dev, &factorized_matrix [0], n * n * sizeof (datatype), cudaMemcpyHostToDevice));
			HANDLE_ERROR (cudaMemcpy (ipiv_dev, &ipiv [0], n * sizeof (int), cudaMemcpyHostToDevice));
		}
		
		template <class datatype>
		void solver <datatype>::execute () {
			TRACE ("Setting up solve...");
			
			utils::copy (n, data_in, data_out);
		
			utils::add_scaled (1, alpha_0 * timestep, rhs, data_out);
			utils::add_scaled (1, alpha_n * timestep, rhs + n - 1, data_out + n - 1);
			utils::add_scaled (n - 2, timestep, rhs + 1, data_out + 1);

			TRACE ("Solving...");

			utils::matrix_solve (n, factorized_matrix_dev, ipiv_dev, data_out);
			
			TRACE ("Solved.");
		}
		
		template class solver <double>;
		template class solver <float>;
	} /* one_d */
} /* cuda */
