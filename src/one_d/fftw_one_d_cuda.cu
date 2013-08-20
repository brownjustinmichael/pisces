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
#include "../utils/utils_cublas.hcu"

namespace one_d
{
	namespace cuda
	{
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
		
		fftw_cosine::fftw_cosine (bases::element* i_element_ptr, int i_n, int i_name_in, int i_name_out) : 
		bases::transform (i_element_ptr, i_n, i_name_in, i_name_out) {
			TRACE ("Instantiating...");
			if (cudaMalloc ((void **) &data_real, 2 * n * sizeof (double)) != cudaSuccess){
				FATAL ("Failed to allocate.\n");
				throw 1;	
			}
			if (cudaMalloc ((void **) &data_complex, n * sizeof (cufftDoubleComplex)) != cudaSuccess) {
				FATAL ("Failed to allocate.\n");
				throw 1;	
			}
			plan = new cufftHandle;
			if (cufftPlan1d(plan, 2 * n - 2, CUFFT_D2Z, 1) != CUFFT_SUCCESS){
				FATAL ("Plan creation failed");
				throw 1;	
			}
			TRACE ("Instantiated.");
			temp.resize (2 * n - 2, 1.4);
			// fourier_plan = fftw_plan_r2r_1d (n, data_in, data_out, FFTW_REDFT00, FFTW_ESTIMATE);
		}
		
		fftw_cosine::~fftw_cosine () {
			cufftDestroy (*plan);
			delete plan;
			cudaFree (data_real);
			cudaFree (data_complex);
		}
		
		void fftw_cosine::execute () {
			bases::transform::execute ();
			if (*flags_ptr & transformed) {
				*flags_ptr &= ~transformed;
			} else {
				*flags_ptr |= transformed;
			}
			
			for (int i = 0; i < n; ++i) {
				DEBUG ("THERE: " << data_in [i]);
			}

			if (cudaMemcpy (data_real, data_in, n * sizeof (double), cudaMemcpyHostToDevice) != cudaSuccess) {
				FATAL ("FAILURE");
				throw 1;
			}
		
			symmetrize <<<1, std::min (n, 512)>>> (n, data_real);
		
			if (cudaThreadSynchronize() != cudaSuccess){
				FATAL ("Failed to synchronize\n");
				throw 1;	
			}
			/* Use the CUFFT plan to transform the signal in place. */
			if (cufftExecD2Z(*plan, data_real, data_complex) != CUFFT_SUCCESS){
				FATAL ("ExecD2Z Forward failed");
				throw 1;	
			}
		
			if (cudaThreadSynchronize() != cudaSuccess){
				FATAL ("Failed to synchronize\n");
				throw 1;	
			}
		
			complex_to_real <<<1, std::min (n, 512)>>> (n, data_complex, data_real);
		
			cudaMemcpy (data_out, data_real, n * sizeof (double), cudaMemcpyDeviceToHost);
		
			if (cudaThreadSynchronize() != cudaSuccess){
				FATAL ("Failed to synchronize\n");
				throw 1;	
			}
		
			for (int i = 0; i < n; ++i) {
				DEBUG ("REAL: " << data_out [i]);
			}
	
			if (cudaThreadSynchronize() != cudaSuccess){
				FATAL ("Failed to synchronize.\n");
				throw 1;	
			}
			
			for (int i = 0; i < n; ++i) {
				data_out [i] *= scalar;
			}
		}
	} /* cuda */
} /* one_d */
