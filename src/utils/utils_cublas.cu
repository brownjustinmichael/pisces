/*!**********************************************************************
 * \file utils_cublas.cu
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "utils_cublas.hcu"
#include <vector>
#include <cassert>

namespace utils
{
	namespace cuda
	{
		void HANDLE_STATUS (cublasStatus status) {
			switch (status) {
				case CUBLAS_STATUS_NOT_INITIALIZED: printf ("CUBLAS didn't initialize correctly.\n"); throw 0; 
				case CUBLAS_STATUS_ALLOC_FAILED: printf ("CUBLAS allocation failed.\n"); throw 0; 
				case CUBLAS_STATUS_INVALID_VALUE: printf ("CUBLAS unsupported value or parameter.\n"); throw 0; 
				case CUBLAS_STATUS_ARCH_MISMATCH: printf ("CUBLAS feature absent in current architecture.\n"); throw 0; 
				case CUBLAS_STATUS_MAPPING_ERROR: printf ("CUBLAS access to GPU memory failed.\n"); throw 0;
				case CUBLAS_STATUS_EXECUTION_FAILED: printf ("CUBLAS failed to execute.\n"); throw 0;
				case CUBLAS_STATUS_INTERNAL_ERROR: printf ("CUBLAS internal operation failed.\n"); throw 0;
			}
		}
		
		config config_instance;
	
		void copy (int n, vect_float* x, vect_float* y, int incx, int incy) {
			cublasScopy (n, x->pointer (), incx, y->pointer (), incy);
		}
	
		void copy (int n, vect_double* x, vect_double* y, int incx, int incy) {
			cublasDcopy (n, x->pointer (), incx, y->pointer (), incy);
		}
	
		void scale (int n, float a, vect_float* x, int incx) {
			cublasSscal (n, a, x->pointer (), incx);
		}
		
		void scale (int n, double a, vect_double* x, int incx) {
			cublasDscal (n, a, x->pointer (), incx);
		}

		double dot (int n, vect_float* x, vect_float* y, int incx, int incy) {
			return cublasSdot (n, x->pointer (), incx, y->pointer (), incy);
		}
	
		double dot (int n, vect_double* x, vect_double* y, int incx, int incy) {
			return cublasDdot (n, x->pointer (), incx, y->pointer (), incy);
		}

		void add_scaled (int n, float a, vect_float* x, vect_float* y, int incx, int incy) {
			cublasSaxpy (n, a, x->pointer (), incx, y->pointer (), incy);
		}
	
		void add_scaled (int n, double a, vect_double* x, vect_double* y, int incx, int incy) {
			cublasDaxpy (n, a, x->pointer (), incx, y->pointer (), incy);
		}
	
		void matrix_vector_multiply (int m, int n, float alpha, vect_float* a, vect_float* x, float beta, vect_float* y, int lda, int incx, int incy) {
			char charN = 'N';
		
			assert (x != y);
		
			if (lda == -1) {
				lda = m;
			}
			
			cublasSgemv (charN, m, n, alpha, a->pointer (), lda, x->pointer (), incx, beta, y->pointer (), incy);
		}
		
		void matrix_vector_multiply (int m, int n, double alpha, vect_double* a, vect_double* x, double beta, vect_double* y, int lda, int incx, int incy) {
			char charN = 'N';
		
			assert (x != y);
		
			if (lda == -1) {
				lda = m;
			}
			
			cublasDgemv (charN, m, n, alpha, a->pointer (), lda, x->pointer (), incx, beta, y->pointer (), incy);
		}
	} /* cublas */
	
} /* utils */
