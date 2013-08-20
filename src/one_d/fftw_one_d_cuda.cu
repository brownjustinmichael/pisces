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
#include "../utils/utils_cublas.cuh"

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
		
		fftw_cosine::fftw_cosine (bases::element* i_element_ptr, int i_n, int i_name_in, int i_name_out) : 
		bases::explicit_plan (i_element_ptr, i_n, i_name_in, i_name_out) {
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
			// fourier_plan = fftw_plan_r2r_1d (n, data_in, data_out, FFTW_REDFT00, FFTW_ESTIMATE);
		}
		
		fftw_cosine::~fftw_cosine () {
			cufftDestroy (*plan);
			delete plan;
			cudaFree (data_real);
			cudaFree (data_complex);
		}
		
		void fftw_cosine::execute () {
			bases::explicit_plan::execute ();
			std::vector <cufftDoubleComplex> temp (n);
			std::vector <double> reverse (2 * n);
			
			for (int i = 0; i < n; ++i) {
				reverse [i] = (double) i;
				reverse [2 * n - i - 2] = (double) i;
			}
			for (int i = 0; i < 2 * n - 2; ++i) {
				DEBUG ("HERE: " << reverse [i]);
			}

			// if (cudaMemcpy (data_real, data_in, n * sizeof (double), cudaMemcpyHostToDevice) != cudaSuccess) {
			// 	FATAL ("FAILURE");
			// 	throw 1;
			// }
			
			if (cudaMemcpy (data_real, &reverse [0], 2 * (n - 1) * sizeof (double), cudaMemcpyHostToDevice) != cudaSuccess) {
				FATAL ("FAILURE");
				throw 1;
			}
		
			// symmetrize <<<1, std::min (n, 512)>>> (n, data_real);
		
			if (cudaThreadSynchronize() != cudaSuccess){
				FATAL ("Failed to synchronize\n");
				throw 1;	
			}
			/* Use the CUFFT plan to transform the signal in place. */
			if (cufftExecD2Z(*plan, data_real, data_complex) != CUFFT_SUCCESS){
				FATAL ("ExecD2Z Forward failed");
				throw 1;	
			}

			cudaMemcpy (&temp [0], data_complex, n * sizeof (cufftDoubleComplex), cudaMemcpyDeviceToHost);
			
			if (cudaThreadSynchronize() != cudaSuccess){
				FATAL ("Failed to synchronize\n");
				throw 1;	
			}
			
			for (int i = 0; i < n; ++i) {
				DEBUG ("REAL: " << temp [i].x);
				DEBUG ("IMAG: " << temp [i].y);
			}		
			
			complex_to_real <<<1, std::min (n, 512)>>> (n, data_complex, data_real);
		
			if (cudaDeviceSynchronize() != cudaSuccess){
				FATAL ("Failed to synchronize\n");
				throw 1;	
			}

			cudaMemcpy (&reverse [0], data_real, n * sizeof (double), cudaMemcpyDeviceToHost);
		
			if (cudaThreadSynchronize() != cudaSuccess){
				FATAL ("Failed to synchronize\n");
				throw 1;	
			}
		
			for (int i = 0; i < n; ++i) {
				DEBUG ("REAL: " << reverse [i]);
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
