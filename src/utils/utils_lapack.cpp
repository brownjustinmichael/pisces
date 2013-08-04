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

namespace utils
{
	void copy (int n, double* dx, double* dy, int incx, int incy) {
		dcopy_ (&n, dx, &incx, dy, &incy);
	}
	
	void scale (int n, double da, double* dx, int incx) {
		dscal_ (&n, &da, dx, &incx);
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
	
	void matrix_factorize (int m, int n, double* a, int *ipiv, int *info, int lda) {
		if (!info) {
			int iinfo;
			info = &iinfo;
		}
		
		if (lda == -1) {
			lda = m;
		}
		
		dgetrf_ (&m, &n, a, &lda, ipiv, info);
	}
	
	void matrix_solve (int n, double* a, int* ipiv, double* b, int *info, int nrhs, int lda, int ldb) {
		char charN = 'N';
		
		if (!info) {
			int iinfo;
			info = &iinfo;
		}
		
		if (lda == -1) {
			lda = n;
		}
		
		if (ldb == -1) {
			ldb = n;
		}
		
		dgetrs_ (&charN, &n, &nrhs, a, &lda, ipiv, b, &ldb, info);
	}
	
	double interpolate (int n, double* dx, double* dy, double x) {
		int i = 1;
		/*
			TODO Allow for reverse dx as well
		*/
		if (x < dx [0] || x > dx [n - 1]) {
			throw 0;
			/*
				TODO better exception?
			*/
		}
		while (x < dx [i]) {
			++i;
		}
		if (x == dx [i]) {
			return dy [i];
		} else {
			return (dy [i] - dy [i - 1]) / (dx [i] - dx [i - 1]) * (x - dx [i]) + dy [i];
		}
	}
	
	void matrix_interpolate (int n, double* dx, int m, double* dy, int incy, double* douty, double x) {
		int i = 1;
		/*
			TODO Allow for reverse dx as well
		*/
		if (x < dx [0] || x > dx [n - 1]) {
			throw 0;
			/*
				TODO better exception?
			*/
		}
		while (x < dx [i]) {
			++i;
		}
		if (x == dx [i]) {
			for (int j = 0; j < m; ++j) {
				douty [j * incy] = dy [j * m + i];
			}
		} else {
			for (int j = 0; j < m; ++j) {
				douty [j * incy] = (dy [j * m + i] - dy [j * m + i - 1]) / (dx [i] - dx [i - 1]) * (x - dx [i]) + dy [j * m + i];
			}
		}
	}
	
} /* utils */