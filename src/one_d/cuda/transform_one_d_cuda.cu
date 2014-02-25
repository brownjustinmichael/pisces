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
#include "../../utils/cuda/utils_cuda.cuh"

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

__global__ void symmetrize (int n, double* data) {
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < n - 1 && tid != 0) {
		data [2 * n - 2 - tid] = data [tid];
		tid += blockDim.x * gridDim.x;
	}
}

__global__ void symmetrize (int n, float* data) {
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < n - 1 && tid != 0) {
		data [2 * n - 2 - tid] = data [tid];
		tid += blockDim.x * gridDim.x;
	}
}

namespace cuda
{
	namespace one_d
	{
		template <>
		transform <double>::transform (int i_n, double* i_data_dev) : 
		n (i_n),
		data_real (i_data_dev) {
			TRACE ("Instantiating...");
			HANDLE_ERROR (cudaMalloc ((void **) &data_complex, n * sizeof (cufftDoubleComplex)));
			cu_plan = new cufftHandle;
			HANDLE_CUFFT (cufftPlan1d((cufftHandle*) cu_plan, 2 * n - 2, CUFFT_D2Z, 1));
			scalar = sqrt (1.0 / 2.0 / ((double) n - 1.0));
			TRACE ("Instantiated.");
		}
		
		template <>
		transform <float>::transform (int i_n, float* i_data_dev) : 
		n (i_n),
		data_real (i_data_dev) {
			TRACE ("Instantiating...");
			HANDLE_ERROR (cudaMalloc ((void **) &data_complex, n * sizeof (cufftComplex)));
			cu_plan = new cufftHandle;
			HANDLE_CUFFT (cufftPlan1d((cufftHandle*) cu_plan, 2 * n - 2, CUFFT_R2C, 1));
			scalar = sqrt (1.0 / 2.0 / ((float) n - 1.0));
			
			HANDLE_ERROR (cudaDeviceSynchronize ());
			
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
			symmetrize <<<1, std::min (n, 512)>>> (n, (double *) data_real);
			
			HANDLE_CUFFT (cufftExecD2Z(*((cufftHandle *) cu_plan), (double*) data_real, (cufftDoubleComplex*) data_complex));
			
			complex_to_real <<<1, std::min (n, 512)>>> (n, (cufftDoubleComplex*) data_complex, (double*) data_real);
			
			utils::scale (n, scalar, (double*) data_real);
		}
		
		template <>
		void transform <float>::execute () {
			symmetrize <<<1, std::min (n, 512)>>> (n, (float *) data_real);
			
			HANDLE_CUFFT (cufftExecR2C(*((cufftHandle *) cu_plan), (float*) data_real, (cufftComplex*) data_complex));
			
			complex_to_real <<<1, std::min (n, 512)>>> (n, (cufftComplex*) data_complex, (float*) data_real);
			
			utils::scale (n, scalar, (float*) data_real);
		}
		
		template <class datatype>
		transfer <datatype>::transfer (int i_n, datatype* i_data_dev, datatype* i_data) :
		n (i_n),
		data_dev (i_data_dev),
		data (i_data) {}
		
		template <class datatype>
		void transfer <datatype>::execute () {
			cudaMemcpy (data, data_dev, n * sizeof (datatype), cudaMemcpyDeviceToHost);
		}
		
		template class transfer <double>;
		template class transfer <float>;
	} /* one_d */
} /* cuda */
