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

#define CUBLAS

namespace utils
{
	template <class datatype>
	struct cublas_vector 
	{
	public:
		cublas_vector (int i_n, datatype *x = NULL, int incx = 1) {
			cublasStatus status;
			inc_vect = 1;
			n = i_n;
			status = cublasAlloc (n, sizeof (datatype), (void**) &vect);
			if (status != CUBLAS_STATUS_SUCCESS) {
				printf ("Couldn't allocate.\n");
				throw 0;
			}
			if (x) {
				copy_to_device (x, incx);
			}
		}
		
		~cublas_vector () {
			cublasStatus status;
			status = cublasFree (vect);
			if (status != CUBLAS_STATUS_SUCCESS) {
				printf ("Couldn't free.\n");
				throw 0;
			}
		}
		
		void copy_to_device (datatype* x, int incx = 1) {
			cublasStatus status;
			status = cublasSetVector (n, sizeof (datatype), x, incx, vect, inc_vect);
			if (status != CUBLAS_STATUS_SUCCESS) {
				printf ("Couldn't copy to device.\n");
				throw 0;
			}
		}
		
		void copy_to_host (datatype* x, int incx = 1) {
			cublasStatus status;
			status = cublasGetVector (n, sizeof (datatype), vect, inc_vect, x, incx);
			if (status != CUBLAS_STATUS_SUCCESS) {
				printf ("Couldn't copy to host.\n");
				throw 0;
			}
		}
		
		datatype* get_pointer () {
			return vect;
		}
		
		datatype* vect;
		int inc_vect;
		int n;
	};
	
	/*
		TODO Add initialize class or function that keeps track of whether the library has initialized and shuts it down properly
		* Also, potentially put all device vectors in there so that they need only be allocated once and are freed at the end
	*/
	
	template <class datatype>
	void copy (int n, datatype* x, datatype* y, int incx, int incy) {
		cublasStatus status;
		
		cublasInit ();
		
		cublas_vector <datatype> vect (n, x, incx);
		
		vect.copy_to_host (y, incy);
		
		status = cublasShutdown();
		if (status != CUBLAS_STATUS_SUCCESS) {
			printf ("Couldn't shut down.");
			throw 0;
		}
	} 
	
	void copy (int n, float* x, float* y, int incx, int incy) {
		copy <float> (n, x, y, incx, incy);
	}
	
	void copy (int n, double* x, double* y, int incx, int incy) {
		copy <double> (n, x, y, incx, incy);
	}
	
	void scale (int n, float a, float* x, int incx) {
		cublasStatus status;
		
		cublasInit ();
		
		cublas_vector <float> vect (n, x, incx);
		
		cublasSscal (n, a, vect.get_pointer (), incx);
		
		vect.copy_to_host (x, incx);
		
		status = cublasShutdown ();
		if (status != CUBLAS_STATUS_SUCCESS) {
			printf ("Couldn't shut down.");
			throw 0;
		}
	}
		
	void scale (int n, double a, double* x, int incx) {
		cublasStatus status;
		
		cublasInit ();
		
		cublas_vector <double> vect (n, x, incx);
		
		cublasDscal (n, a, vect.get_pointer (), incx);
		
		vect.copy_to_host (x, incx);
		
		status = cublasShutdown ();
		if (status != CUBLAS_STATUS_SUCCESS) {
			printf ("Couldn't shut down.");
			throw 0;
		}
	}
	
	double dot (int n, double* dx, double* dy, int incx, int incy) {return 0.0;}
	
	void add_scaled (int n, double da, double *dx, double *dy, int incx, int incy) {}
	
	void matrix_vector_multiply (int m, int n, double alpha, double *a, double *x, double beta, double *y, int lda, int incx, int incy) {}
	
} /* utils */
