/*!**********************************************************************
 * \file utils_cublas.cu
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "utils.hpp"
#include "utils_cublas.hcu"
#include <vector>

namespace utils
{
	namespace cublas
	{
		cublas_config cublas_config_instance;
	
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
	} /* cublas */
	
} /* utils */
