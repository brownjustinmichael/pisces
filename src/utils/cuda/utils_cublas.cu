/*!**********************************************************************
 * \file utils_cublas.cu
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "../../config.hpp"
#include "utils_cublas.hpp"
#include "utils_cuda.hpp"
#include <vector>
#include <cassert>
#include <cublas.h>
#include <stdio.h>

namespace utils
{
	namespace cuda
	{
		namespace cublas
		{
			struct config
			{
			public:
				config () {
					CUBLAS_HANDLE_ERROR (cublasInit ());
				}
	
				virtual ~config () {
					CUBLAS_HANDLE_ERROR (cublasShutdown ());
				}
			};
	
			void copy (int n, float* x, float* y, int incx, int incy) {
				cublasScopy (n, x, incx, y, incy);
			}

			void copy (int n, double* x, double* y, int incx, int incy) {
				cublasDcopy (n, x, incx, y, incy);
			}

			void scale (int n, float a, float* x, int incx) {
				cublasSscal (n, a, x, incx);
			}
	
			void scale (int n, double a, double* x, int incx) {
				cublasDscal (n, a, x, incx);
			}

			double dot (int n, float* x, float* y, int incx, int incy) {
				return cublasSdot (n, x, incx, y, incy);
			}

			double dot (int n, double* x, double* y, int incx, int incy) {
				return cublasDdot (n, x, incx, y, incy);
			}

			void add_scaled (int n, float a, float* x, float* y, int incx, int incy) {
				cublasSaxpy (n, a, x, incx, y, incy);
			}

			void add_scaled (int n, double a, double* x, double* y, int incx, int incy) {
				cublasDaxpy (n, a, x, incx, y, incy);
			}

			void matrix_vector_multiply (int m, int n, float alpha, float* a, float* x, float beta, float* y, int lda, int incx, int incy) {
				char charN = 'N';
	
				assert (x != y);
	
				if (lda == -1) {
					lda = m;
				}
		
				cublasSgemv (charN, m, n, alpha, a, lda, x, incx, beta, y, incy);
			}
	
			void matrix_vector_multiply (int m, int n, double alpha, double* a, double* x, double beta, double* y, int lda, int incx, int incy) {
				char charN = 'N';
	
				assert (x != y);
	
				if (lda == -1) {
					lda = m;
				}
		
				cublasDgemv (charN, m, n, alpha, a, lda, x, incx, beta, y, incy);
			}
	
			config config_instance;
		} /* cublas */
	} /* cuda */
} /* utils */

