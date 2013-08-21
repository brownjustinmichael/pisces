/*!**********************************************************************
 * \file fftw_one_d_cuda.cu
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-16.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <math.h>
#include <cufft.h>
#include "fftw_one_d_cuda.hpp"
#include "../../utils/utils_cublas.hpp"

#define HANDLE_ERROR(status) \
{cudaError_t result = status; \
switch (result) { \
	case cudaErrorMemoryAllocation: FATAL ("Memory Allocation Error."); throw 0; \
	case cudaErrorInvalidValue: FATAL ("Invalid value passed."); throw 0; \
	default: if (status != cudaSuccess) {FATAL ("Other problem."); throw 0;}}}

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

__global__ void symmetrize (int n, double* data) {
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < n - 1 && tid != 0) {
		data [2 * n - 2 - tid] = data [tid];
		tid += blockDim.x * gridDim.x;
	}
}

namespace one_d
{
	namespace cuda
	{
		fftw_cosine::fftw_cosine (bases::element <double>* i_element_ptr, int i_n, double* i_data_dev) : 
		bases::plan <double> (i_element_ptr),
		n (i_n),
		data_real (i_data_dev) {
			TRACE ("Instantiating...");
			HANDLE_ERROR (cudaMalloc ((void **) &data_complex, n * sizeof (cufftDoubleComplex)));
			cu_plan = new cufftHandle;
			HANDLE_CUFFT (cufftPlan1d((cufftHandle*) cu_plan, 2 * n - 2, CUFFT_D2Z, 1));
			scalar = sqrt (1.0 / 2.0 / ((double) n - 1.0));
			TRACE ("Instantiated.");
		}
		
		fftw_cosine::~fftw_cosine () {
			cufftDestroy (*((cufftHandle *) cu_plan));
			delete (cufftHandle*) cu_plan;
			cudaFree (data_complex);
		}
		
		void fftw_cosine::execute () {
			bases::plan <double>::execute ();
			
			symmetrize <<<1, std::min (n, 512)>>> (n, (double *) data_real);
			
			HANDLE_CUFFT (cufftExecD2Z(*((cufftHandle *) cu_plan), (double*) data_real, (cufftDoubleComplex*) data_complex));
			
			complex_to_real <<<1, std::min (n, 512)>>> (n, (cufftDoubleComplex*) data_complex, (double*) data_real);
			
			utils::cuda::scale (n, scalar, (double*) data_real);
		}
		
		transfer::transfer (bases::element <double>* i_element_ptr, int i_n, double* i_data_dev, double* i_data) :
		bases::plan <double> (i_element_ptr),
		data_dev (i_data_dev),
		data (i_data) {}
		
		void transfer::execute () {
			cudaMemcpy (data, data_dev, n * sizeof (double), cudaMemcpyDeviceToHost);
		}
	} /* cuda */
} /* one_d */
