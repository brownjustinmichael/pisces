/*!***********************************************************************
 * \file utils_lapack.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-29.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cassert>
#include "utils.hpp"
#include "lapack.hpp"

#define LAPACK

namespace utils
{
	void copy (int n, float* x, float* y, int incx, int incy) {
		scopy_ (&n, x, &incx, y, &incy);
	}
	
	void copy (int n, double* x, double* y, int incx, int incy) {
		dcopy_ (&n, x, &incx, y, &incy);
	}
	
	void scale (int n, float a, float* x, int incx) {
		sscal_ (&n, &a, x, &incx);
	}
	
	void scale (int n, double a, double* x, int incx) {
		dscal_ (&n, &a, x, &incx);
	}
	
	double dot (int n, double* dx, double* dy, int incx, int incy) {
		return ddot_ (&n, dx, &incx, dy, &incy);
	}	
	
	void add_scaled (int n, double da, double *dx, double *dy, int incx, int incy) {
		daxpy_ (&n, &da, dx, &incx, dy, &incy);
	}
	
	void matrix_vector_multiply (int m, int n, double alpha, double *a, double *x, double beta, double *y, int lda, int incx, int incy) {
		char charN = 'N';
		
		assert (x != y);
		
		if (lda == -1) {
			lda = m;
		}
				
		dgemv_ (&charN, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	
} /* utils */