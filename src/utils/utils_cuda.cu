/*!**********************************************************************
 * \file utils_cuda.cu
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "utils.hpp"
#include <cublas.h>
#include <stdio.h>
#include <vector>

namespace utils
{
	void HANDLE_STATUS (cublasStatus status) {
		switch (status) {
			case CUBLAS_STATUS_NOT_INITIALIZED: printf ("CUBLAS didn't initialize correctly."); throw 0; 
			case CUBLAS_STATUS_ALLOC_FAILED: printf ("CUBLAS allocation failed."); throw 0; 
			case CUBLAS_STATUS_INVALID_VALUE: printf ("CUBLAS unsupported value or parameter."); throw 0; 
			case CUBLAS_STATUS_ARCH_MISMATCH: printf ("CUBLAS feature absent in current architecture."); throw 0; 
			case CUBLAS_STATUS_MAPPING_ERROR: printf ("CUBLAS access to GPU memory failed."); throw 0;
			case CUBLAS_STATUS_EXECUTION_FAILED: printf ("CUBLAS failed to execute."); throw 0;
			case CUBLAS_STATUS_INTERNAL_ERROR: printf ("CUBLAS internal operation failed."); throw 0;
		}
	}
	template <class datatype>
	struct cublas_vector 
	{
	public:
		cublas_vector (int i_n, datatype *x = NULL, int incx = 1) {
			inc_vect = 1;
			n = i_n;
			is_active = true;
			HANDLE_STATUS (cublasAlloc (n, sizeof (datatype), (void**) &vect));
			if (x) {
				copy_to_device (n, x, incx);
			}
		}
		
		~cublas_vector () {
			HANDLE_STATUS (cublasFree (vect));
		}
		
		datatype* operator& () {
			return vect;
		}
		
		int size () {
			return n;
		}
		
		void copy_to_device (int n, datatype* x, int incx = 1) {
			HANDLE_STATUS (cublasSetVector (n, sizeof (datatype), x, incx, vect, inc_vect));
		}
		
		void copy_to_host (int n, datatype* x, int incx = 1) {
			HANDLE_STATUS (cublasGetVector (n, sizeof (datatype), vect, inc_vect, x, incx));
		}
		
	private:
		int n;
		datatype* vect;
		int inc_vect;
	};
	
	
	struct cublas_handler
	{
	public:
		cublas_handler () {
			HANDLE_STATUS (cublasInit ());
		}
		
		virtual ~cublas_handler () {
			HANDLE_STATUS (cublasShutdown ());
		}
	};
	
	cublas_handler handler;
	
	/*
		TODO Also, potentially put all device vectors in handler so that they need only be allocated once and are freed at the end
	*/
	
	void copy (int n, float* x, float* y, int incx, int incy) {
		cublas_vector <float> vect (n, x, incx);
		vect.copy_to_host (n, y, incy);
	}
	
	void copy (int n, double* x, double* y, int incx, int incy) {
		cublas_vector <double> vect (n, x, incx);
		vect.copy_to_host (n, y, incy);
	}
	
	void scale (int n, float a, float* x, int incx) {
		cublas_vector <float> vect (n, x, incx);
		cublasSscal (n, a, &vect, incx);
		vect.copy_to_host (n, x, incx);
	}
		
	void scale (int n, double a, double* x, int incx) {
		cublas_vector <double> vect (n, x, incx);
		cublasDscal (n, a, &vect, incx);
		vect.copy_to_host (n, x, incx);
	}
	
	double dot (int n, double* dx, double* dy, int incx, int incy) {return 0.0;}
	
	void add_scaled (int n, double da, double *dx, double *dy, int incx, int incy) {}
	
	void matrix_vector_multiply (int m, int n, double alpha, double *a, double *x, double beta, double *y, int lda, int incx, int incy) {}
	
} /* utils */
