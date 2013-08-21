/*!**********************************************************************
 * \file utils_cublas.cu
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#define HANDLE_STATUS(status) \
{cublasStatus result = status; \
switch (result) { \
	case CUBLAS_STATUS_NOT_INITIALIZED: printf ("CUBLAS didn't initialize correctly.\n"); throw 0; \
	case CUBLAS_STATUS_ALLOC_FAILED: printf ("CUBLAS allocation failed.\n"); throw 0; \
	case CUBLAS_STATUS_INVALID_VALUE: printf ("CUBLAS unsupported value or parameter.\n"); throw 0; \
	case CUBLAS_STATUS_ARCH_MISMATCH: printf ("CUBLAS feature absent in current architecture.\n"); throw 0; \
	case CUBLAS_STATUS_MAPPING_ERROR: printf ("CUBLAS access to GPU memory failed.\n"); throw 0; \
	case CUBLAS_STATUS_EXECUTION_FAILED: printf ("CUBLAS failed to execute.\n"); throw 0; \
	case CUBLAS_STATUS_INTERNAL_ERROR: printf ("CUBLAS internal operation failed.\n"); throw 0;}}

#include "utils_cublas.hpp"
#include <vector>
#include <cassert>
#include <cublas.h>
#include <stdio.h>

namespace utils
{
	namespace cuda
	{
		template <class datatype>
		vector <datatype>::vector (int i_n, datatype *x, int incx) {
			inc_vect = 1;
			n = i_n;
			if (n != 0) {
				HANDLE_STATUS (cublasAlloc (n, sizeof (datatype), (void**) &vect));
				if (x) {
					copy_to_device (n, x, incx);
				}
			}
		}
		
		template <class datatype>
		vector <datatype>::~vector () {
			if (n != 0) {
				HANDLE_STATUS (cublasFree (vect));
			}
		}
		
		template <class datatype>
		void vector <datatype>::resize (int i_n) {
			if (n != 0) {
				HANDLE_STATUS (cublasFree (vect));
			}
			if (i_n != 0) {
				n = i_n;
				HANDLE_STATUS (cublasAlloc (n, sizeof (datatype), (void**) &vect));
			}
		}
		
		template <class datatype>
		void vector <datatype>::copy_to_device (int i_n, datatype* x, int incx) {
			if (i_n != 0 && n != 0) {
				HANDLE_STATUS (cublasSetVector (i_n, sizeof (datatype), x, incx, vect, inc_vect));
			} else if (n == 0) {
				throw 0;
			}
		}
	
		template <class datatype>
		void vector <datatype>::copy_to_host (int i_n, datatype* x, int incx) {
			if (i_n != 0 && n != 0) {
				HANDLE_STATUS (cublasGetVector (i_n, sizeof (datatype), vect, inc_vect, x, incx));
			} else if (n == 0) {
				throw 0;
			}
		}
		
		struct config
		{
		public:
			config () {
				HANDLE_STATUS (cublasInit ());
			}
		
			virtual ~config () {
				HANDLE_STATUS (cublasShutdown ());
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
		
		template class vector <double>;
		template class vector <float>;
	} /* cuda */
} /* utils */
