/*!**********************************************************************
 * \file transform_two_d.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-15.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cmath>
#include "transform_two_d_cuda.hpp"

namespace two_d
{
	namespace fourier
	{
		template <>
		horizontal_transform <float>::horizontal_transform (int n, int m, float* i_data_in, float* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags) :
		plans::plan <float> (i_element_flags, i_component_flags),
		n (n), 
		m (m), 
		data_in (i_data_in),
		data_out (i_data_out ? i_data_out : i_data_in),
		flags (i_flags) {
			TRACE ("Initializing...");
			scalar = 1.0 / std::sqrt (n);
			cu_plan = new cufftHandle;
			
			if (!(flags & inverse)) {
				cufftPlanMany ((cufftHandle *) cu_plan, &m, &ldm, ldn, 1, &ldm, ldn, 1, CFFT_R2C, n);
			} else {
				cufftPlanMany ((cufftHandle *) cu_plan, &m, &ldm, ldn, 1, &ldm, ldn, 1, CFFT_C2R, n);
			}
		}
		
		template <>
		horizontal_transform <float>::horizontal_transform (plans::grid <float> &i_grid_n, plans::grid <float> &i_grid_m, float* i_data_in, float* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags) :
		plans::plan <float> (i_element_flags, i_component_flags),
		n (i_grid_n.n), 
		m (i_grid_m.n), 
		data_in (i_data_in),
		data_out (i_data_out ? i_data_out : i_data_in),
		flags (i_flags) {
			TRACE ("Initializing...");
			scalar = 1.0 / std::sqrt (n);
			cu_plan = new cufftHandle;
			
			if (!(flags & inverse)) {
				cufftPlanMany ((cufftHandle *) cu_plan, &m, &ldm, ldn, 1, &ldm, ldn, 1, CFFT_R2C, n);
			} else {
				cufftPlanMany ((cufftHandle *) cu_plan, &m, &ldm, ldn, 1, &ldm, ldn, 1, CFFT_R2C, n);
			}
		}
		
		template <>
		horizontal_transform <float>::~horizontal_transform () {
			cufftDestroy (cu_plan);
		}
		
		template <>
		horizontal_transform <double>::horizontal_transform (int n, int m, double* i_data_in, double* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags) :
		plans::plan <double> (i_element_flags, i_component_flags),
		n (n), 
		m (m), 
		data_in (i_data_in),
		data_out (i_data_out ? i_data_out : i_data_in),
		flags (i_flags) {
			TRACE ("Initializing...");
			scalar = 1.0 / std::sqrt (n);
			cu_plan = new cufftHandle;
			
			if (!(flags & inverse)) {
				cufftPlanMany ((cufftHandle *) cu_plan, &m, &ldm, ldn, 1, &ldm, ldn, 1, CFFT_D2Z, n);
			} else {
				cufftPlanMany ((cufftHandle *) cu_plan, &m, &ldm, ldn, 1, &ldm, ldn, 1, CFFT_Z2D, n);
			}
		}
		
		template <>
		horizontal_transform <double>::horizontal_transform (plans::grid <double> &i_grid_n, plans::grid <double> &i_grid_m, double* i_data_in, double* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags) :
		plans::plan <double> (i_element_flags, i_component_flags),
		n (i_grid_n.n), 
		m (i_grid_m.n), 
		data_in (i_data_in),
		data_out (i_data_out ? i_data_out : i_data_in),
		flags (i_flags) {
			TRACE ("Initializing...");
			scalar = 1.0 / std::sqrt (n);
			cu_plan = new cufftHandle;
			
			if (!(flags & inverse)) {
				cufftPlanMany ((cufftHandle *) cu_plan, &m, &ldm, ldn, 1, &ldm, ldn, 1, CFFT_D2Z, n);
			} else {
				cufftPlanMany ((cufftHandle *) cu_plan, &m, &ldm, ldn, 1, &ldm, ldn, 1, CFFT_Z2D, n);
			}
		}
		
		template <>
		horizontal_transform <double>::~horizontal_transform () {
			cufftDestroy (cu_plan);
		}
		
		template <>
		void horizontal_transform <float>::execute () {
			TRACE ("Executing...");
			
			if (!(flags & inverse)) {
				cufftExecR2C ((cufftHandle *) cu_plan, data_in, data_out);
			} else {
				cufftExecC2R ((cufftHandle *) cu_plan, data_in, data_out);
			}
			
			if (*component_flags & transformed_horizontal) {
				*component_flags &= ~transformed_horizontal;
			} else {
				*component_flags |= transformed_horizontal;
			}
				
			utils::cuda::cublas::scale (2 * (n / 2 + 1) * m, scalar, (float*) data_out);
		}
		
		template <>
		void horizontal_transform <double>::execute () {
			TRACE ("Executing...");
			
			if (!(flags & inverse)) {
				cufftExecD2Z ((cufftHandle *) cu_plan, data_in, data_out);
			} else {
				cufftExecZ2D ((cufftHandle *) cu_plan, data_in, data_out);
			}
			
			if (*component_flags & transformed_horizontal) {
				*component_flags &= ~transformed_horizontal;
			} else {
				*component_flags |= transformed_horizontal;
			}
		
			utils::cuda::cublas::scale (2 * (n / 2 + 1) * m, scalar, (double*) data_out);
			
			TRACE ("Execution Complete.");
		}
		
		template class horizontal_transform <float>;
		template class horizontal_transform <double>;
		
		template <>
		vertical_transform <float>::vertical_transform (int n, int m, float* i_data_in, float* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags) :
		plans::plan <float> (i_element_flags, i_component_flags),
		n (n), 
		m (m), 
		data_in (i_data_in),
		data_out (i_data_out ? i_data_out : i_data_in),
		flags (i_flags) {
			TRACE ("Initializing...");
			scalar = 1.0;
			if (m > 1 && !(flags & ignore_m)) {
				scalar /= std::sqrt (2.0 * (m - 1));
			}
			int ldn = 2 * n - 2;
			cu_plan = new cufftHandle;
			HANDLE_ERROR (cudaMalloc ((void **) &data_complex, (n - 1) * m * sizeof (cufftDoubleComplex)));
			HANDLE_ERROR (cudaMalloc ((void **) &data_real, 2 * (n - 1) * m * sizeof (double)));

			cufftPlanMany ((cufftHandle *) cu_plan, 1, &ldn, NULL, 1, ldn, NULL, 1, ldn / 2, CUFFT_R2C, m);
			TRACE ("Instantiated.");
		}
		
		template <>
		vertical_transform <float>::vertical_transform (plans::grid <float> &i_grid_n, plans::grid <float> &i_grid_m, float* i_data_in, float* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags) :
		plans::plan <float> (i_element_flags, i_component_flags),
		n (i_grid_n.n), 
		m (i_grid_m.n), 
		data_in (i_data_in),
		data_out (i_data_out ? i_data_out : i_data_in),
		flags (i_flags) {
			TRACE ("Initializing...");
			scalar = 1.0;
			if (m > 1 && !(flags & ignore_m)) {
				scalar /= std::sqrt (2.0 * (m - 1));
			}
			int ldn = 2 * n - 2;
			cu_plan = new cufftHandle;
			HANDLE_ERROR (cudaMalloc ((void **) &data_complex, (n - 1) * m * sizeof (cufftDoubleComplex)));
			HANDLE_ERROR (cudaMalloc ((void **) &data_real, 2 * (n - 1) * m * sizeof (double)));

			cufftPlanMany ((cufftHandle *) cu_plan, 1, &ldn, NULL, 1, ldn, NULL, 1, ldn / 2, CUFFT_R2C, m);
			TRACE ("Instantiated.");
		}
		
		template <>
		vertical_transform <float>::~vertical_transform () {
			fftwf_destroy_plan (plans_float [0]);
		}
		
		template <>
		vertical_transform <double>::vertical_transform (int n, int m, double* i_data_in, double* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags) :
		plans::plan <double> (i_element_flags, i_component_flags),
		n (n), 
		m (m), 
		data_in (i_data_in),
		data_out (i_data_out ? i_data_out : i_data_in),
		flags (i_flags) {
			TRACE ("Initializing...");
			scalar = 1.0;
			if (m > 1 && !(flags & ignore_m)) {
				scalar /= std::sqrt (2.0 * (m - 1));
			}
			
			fftw_r2r_kind kind = FFTW_REDFT00;

			plans.resize (threads);
			
			int index = 0;
			for (int i = 0; i < threads; ++i) {
				int nn = (2 * (n / 2 + 1)) / threads + (i < (2 * (n / 2 + 1) % threads)? 1: 0);
				plans [i] = fftw_plan_many_r2r (1, &m, nn, data_in + index * m, NULL, 1, m, data_out + index * m, NULL, 1, m, &kind, FFTW_MEASURE);
				index += nn;
			}
		}

		template <>
		vertical_transform <double>::vertical_transform (plans::grid <double> &i_grid_n, plans::grid <double> &i_grid_m, double* i_data_in, double* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags) :
		plans::plan <double> (i_element_flags, i_component_flags),
		n (i_grid_n.n), 
		m (i_grid_m.n), 
		data_in (i_data_in),
		data_out (i_data_out ? i_data_out : i_data_in),
		flags (i_flags) {
			TRACE ("Initializing...");
			scalar = 1.0;
			if (m > 1 && !(flags & ignore_m)) {
				scalar /= std::sqrt (2.0 * (m - 1));
			}
			
			fftw_r2r_kind kind = FFTW_REDFT00;

			plans.resize (threads);
			
			int index = 0;
			for (int i = 0; i < threads; ++i) {
				int nn = (2 * (n / 2 + 1)) / threads + (i < (2 * (n / 2 + 1) % threads)? 1: 0);
				plans [i] = fftw_plan_many_r2r (1, &m, nn, data_in + index * m, NULL, 1, m, data_out + index * m, NULL, 1, m, &kind, FFTW_MEASURE);
				index += nn;
			}
		}
		
		template <>
		vertical_transform <double>::~vertical_transform () {
			fftw_destroy_plan (plans [0]);
		}
		
		
		template <>
		void vertical_transform <float>::execute () {
			TRACE ("Executing...");
	
			if (m > 1 && !(flags & ignore_m)) {
				// #pragma omp parallel for num_threads (threads)
				for (int i = 0; i < threads; ++i) {
					fftwf_execute (plans_float [i]);
				}
			}
			
			if (*component_flags & transformed_vertical) {
				*component_flags &= ~transformed_vertical;
			} else {
				*component_flags |= transformed_vertical;
			}
								
			for (int i = 0; i < 2 * (n / 2 + 1) * m; ++i) {
				data_out [i] *= scalar;
			}
		}
		
		template <>
		void vertical_transform <double>::execute () {
			TRACE ("Executing...");

			if (m > 1 && !(flags & ignore_m)) {
				// #pragma omp parallel for num_threads (threads)
				for (int i = 0; i < threads; ++i) {
					fftw_execute (plans [i]);
				}
			}
			
			if (*component_flags & transformed_vertical) {
				*component_flags &= ~transformed_vertical;
			} else {
				*component_flags |= transformed_vertical;
			}
			
			for (int i = 0; i < 2 * (n / 2 + 1) * m; ++i) {
				data_out [i] *= scalar;
			}
			
			TRACE ("Execution Complete.");
		}
		
		template class vertical_transform <float>;
		template class vertical_transform <double>;
	} /* fourier */
} /* two_d */

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
		transform <double>::transform (plans::grid <double> &i_grid, double* i_data_in, double* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags) : 
		plans::plan <double> (i_element_flags, i_component_flags), 
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
		transform <float>::transform (plans::grid <float> &i_grid, float* i_data_in, float* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags) : 
		plans::plan <float> (i_element_flags, i_component_flags), 
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
