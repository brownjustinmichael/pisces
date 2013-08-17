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
#define PI 3.141592653589793

namespace one_d
{
	namespace cuda
	{
		__global__ void real_to_complex (int n, double* in, cufftComplex* out) {
			int tid = threadIdx.x;
			out [tid].x = in [tid];
			out [tid].y = 0.0;
		}
		
		fftw_cosine::fftw_cosine (bases::element* i_element_ptr, int i_n, int i_name_in, int i_name_out) : 
		bases::transform (i_element_ptr, i_n, i_name_in, i_name_out) {
			cudaMalloc ((void **) data, n * sizeof (cufftComplex));
			if (cudaGetLastError() != cudaSuccess){
				fprintf(stderr, "Cuda error: Failed to allocate\n");
				throw 1;	
			}
			
			// fourier_plan = fftw_plan_r2r_1d (n, data_in, data_out, FFTW_REDFT00, FFTW_ESTIMATE);
		}
		
		fftw_cosine::~fftw_cosine () {
			
		}
		
		void fftw_cosine::execute () {
			bases::transform::execute ();
			if (*flags_ptr & transformed) {
				*flags_ptr &= ~transformed;
			} else {
				*flags_ptr |= transformed;
			}
			
			// fftw_execute (fourier_plan);	

			for (int i = 0; i < n + 1; ++i) {
				data_out [i] *= scalar;
			}
		}
	} /* cuda */
} /* one_d */
