/*!**********************************************************************
 * \file transform_one_d_cuda.cu
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-08-16.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <math.h>
#include <cufft.h>
#include "transform_one_d_cuda.hpp"
#include "../../utils/cuda/utils_cublas.hpp"
#include "../../utils/cuda/utils_cuda.hpp"

#define HANDLE_CUFFT(status) \
{cufftResult result = status; \
switch (result) { \
	case CUFFT_INVALID_PLAN: FATAL ("Invalid plan for cufft."); throw 0; \
	case CUFFT_INVALID_VALUE: FATAL ("Invalid value for cufft."); throw 0; \
	case CUFFT_INTERNAL_ERROR: FATAL ("Internal driver error for cufft."); throw 0; \
	case CUFFT_EXEC_FAILED: FATAL ("Failed to execute transform on cufft."); throw 0; \
	default: if (status != CUFFT_SUCCESS) {FATAL ("Cufft Other problem."); throw 0;}}}

__global__ void real_to_complex (int n, double* in, cufftDoubleComplex* out) {
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < n) {
		out [tid].x = in [tid];
		out [tid].y = 0.0;
		tid += blockDim.x * gridDim.x;
	}
}

__global__ void complex_to_real (int n, cufftDoubleComplex* in, double* out) {
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < n) {
		out [tid] = in [tid].x;
		tid += blockDim.x * gridDim.x;
	}
}

__global__ void complex_to_real (int n, cufftComplex* in, float* out) {
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < n) {
		out [tid] = in [tid].x;
		tid += blockDim.x * gridDim.x;
	}
}

__global__ void symmetrize (int n, double* data_in, double *data_out) {
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < n - 1 && tid != 0) {
		data_out [tid] = data_in [tid];
		data_out [2 * n - 2 - tid] = data_in [tid];
		tid += blockDim.x * gridDim.x;
	}
}

__global__ void symmetrize (int n, float* data_in, float *data_out) {
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < n - 1 && tid != 0) {
		data_out [tid] = data_in [tid];
		data_out [2 * n - 2 - tid] = data_in [tid];
		tid += blockDim.x * gridDim.x;
	}
}

namespace one_d
{
	namespace cuda
	{
		template <>
		transform <double>::transform (bases::grid <double> &i_grid, double* i_data_in, double* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags) : 
		bases::plan <double> (i_element_flags, i_component_flags), 
		n (i_grid.n),
		data_in (i_data_in),
		data_out (i_data_out) {
			TRACE ("Instantiating...");
			HANDLE_ERROR (cudaMalloc ((void **) &data_complex, n * sizeof (cufftDoubleComplex)));
			cu_plan = new cufftHandle;
			HANDLE_CUFFT (cufftPlan1d((cufftHandle*) cu_plan, 2 * n - 2, CUFFT_D2Z, 1));
			scalar = sqrt (1.0 / 2.0 / ((double) n - 1.0));
			TRACE ("Instantiated.");
		}
	
		template <>
		transform <float>::transform (bases::grid <float> &i_grid, float* i_data_in, float* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags) : 
		bases::plan <float> (i_element_flags, i_component_flags), 
		n (i_grid.n),
		data_in (i_data_in),
		data_out (i_data_out) {
			TRACE ("Instantiating...");
			HANDLE_ERROR (cudaMalloc ((void **) &data_complex, n * sizeof (cufftComplex)));
			cu_plan = new cufftHandle;
			HANDLE_CUFFT (cufftPlan1d((cufftHandle*) cu_plan, 2 * n - 2, CUFFT_R2C, 1));
			scalar = sqrt (1.0 / 2.0 / ((float) n - 1.0));
			TRACE ("Instantiated.");
		}
	
		template <class datatype>
		transform <datatype>::~transform () {
			cufftDestroy (*((cufftHandle *) cu_plan));
			delete (cufftHandle*) cu_plan;
			cudaFree (data_complex);
		}
	
		template <>
		void transform <double>::execute () {
			TRACE ("Executing");
			symmetrize <<<1, std::min (n, 512)>>> (n, (double *) data_in, (double *) data_out);
			
			HANDLE_CUFFT (cufftExecD2Z(*((cufftHandle *) cu_plan), (double*) data_out, (cufftDoubleComplex*) data_complex));
			
			complex_to_real <<<1, std::min (n, 512)>>> (n, (cufftDoubleComplex*) data_complex, (double*) data_out);
			
			utils::cuda::cublas::scale (n, scalar, (double*) data_out);
			
			if (!(*component_flags & transformed_vertical)) {
				*component_flags |= transformed_vertical;
			} else {
				*component_flags &= ~transformed_vertical;
			}
		}
	
		template <>
		void transform <float>::execute () {
			TRACE ("Executing");
			symmetrize <<<1, std::min (n, 512)>>> (n, (float *) data_in, (float *) data_out);
			
			HANDLE_CUFFT (cufftExecR2C(*((cufftHandle *) cu_plan), (float*) data_out, (cufftComplex*) data_complex));
			
			complex_to_real <<<1, std::min (n, 512)>>> (n, (cufftComplex*) data_complex, (float*) data_out);
			
			utils::cuda::cublas::scale (n, scalar, (float*) data_out);
			
			if (!(*component_flags & transformed_vertical)) {
				*component_flags |= transformed_vertical;
			} else {
				*component_flags &= ~transformed_vertical;
			}
		}
	} /* cuda */
} /* one_d */
