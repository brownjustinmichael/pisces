/*!***********************************************************************
 * \file utils_lapack.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-29.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cassert>
#include "utils.hpp"

/*!*******************************************************************
 * \brief Function from BLAS that calculates a dot product
 * 
 * \param n A pointer to an integer number of elements in dx and dy to dot
 * \param x The first double array of data
 * \param incx A pointer to an integer spacing of elements in dx
 * \param y The second double array of data
 * \param incy A pointer to an integer spacing of elements in dy
 *********************************************************************/
extern "C" double sdot_ (int *n, float *x, int *incx, float *y, int *incy);

/*!*******************************************************************
 * \brief Function from BLAS that calculates a dot product
 * 
 * \param n A pointer to an integer number of elements in dx and dy to dot
 * \param x The first double array of data
 * \param incx A pointer to an integer spacing of elements in dx
 * \param y The second double array of data
 * \param incy A pointer to an integer spacing of elements in dy
 *********************************************************************/
extern "C" double ddot_ (int *n, double *x, int *incx, double *y, int *incy);

/*!*******************************************************************
 * \brief Function from BLAS that copies a float array to another in place
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param x The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 * \param y The array to which the data are copied
 * \param incy A pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void scopy_ (int *n, float *x, int *incx, float *y, int *incy);

/*!*******************************************************************
 * \brief Function from BLAS that copies a double array to another in place
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param x The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 * \param y The array to which the data are copied
 * \param incy A pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void dcopy_ (int *n, double *x, int *incx, double *y, int *incy);

/*!*******************************************************************
 * \brief Function from BLAS that scales a float array
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param a The float by which to scale the data
 * \param x The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 *********************************************************************/
extern "C" void sscal_ (int *n, float *a, float *x, int *incx);

/*!*******************************************************************
 * \brief Function from BLAS that scales a double array
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param a The double by which to scale the data
 * \param x The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 *********************************************************************/
extern "C" void dscal_ (int *n, double *a, double *x, int *incx);

/*!*******************************************************************
 * \brief Function from BLAS for vector-vector addition (dy = da * dx + dy)
 * 
 * \param n A pointer to the integer number of values in dy
 * \param a A pointer to the double da
 * \param x The double vector dx
 * \param incx A pointer to the integer spacing of elements in dx
 * \param y The double vector dy
 * \param incy A pointer to the integer spacing of elements in dy
 *********************************************************************/
extern "C" void saxpy_ (int *n, float *a, float *x, int *incx, float *y, int *incy);

/*!*******************************************************************
 * \brief Function from BLAS for vector-vector addition (dy = da * dx + dy)
 * 
 * \param n A pointer to the integer number of values in dy
 * \param a A pointer to the double da
 * \param x The double vector dx
 * \param incx A pointer to the integer spacing of elements in dx
 * \param y The double vector dy
 * \param incy A pointer to the integer spacing of elements in dy
 *********************************************************************/
extern "C" void daxpy_ (int *n, double *a, double *x, int *incx, double *y, int *incy);

/*!*******************************************************************
 * \brief Function from BLAS for matrix-vector multiplication (y = alpha * a * x + beta * y)
 * 
 * \param trans A pointer to transposition character ("N" for not transposed, "T" for transposed)
 * \param m A pointer to the number of rows in a
 * \param n A pointer to the number of columns in a
 * \param alpha A pointer to the double multiplier on a
 * \param a The double matrix a
 * \param lda A pointer to the integer number of leading dimension of a
 * \param x The double vector x
 * \param incx A pointer to an integer spacing of elements in x
 * \param beta A pointer to the double multiplier on y
 * \param y The double vector y, overwritten with the solution
 * \param incy A pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void sgemv_ (char *trans, int *m, int *n, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y, int *incy);

/*!*******************************************************************
 * \brief Function from BLAS for matrix-vector multiplication (y = alpha * a * x + beta * y)
 * 
 * \param trans A pointer to transposition character ("N" for not transposed, "T" for transposed)
 * \param m A pointer to the number of rows in a
 * \param n A pointer to the number of columns in a
 * \param alpha A pointer to the double multiplier on a
 * \param a The double matrix a
 * \param lda A pointer to the integer number of leading dimension of a
 * \param x The double vector x
 * \param incx A pointer to an integer spacing of elements in x
 * \param beta A pointer to the double multiplier on y
 * \param y The double vector y, overwritten with the solution
 * \param incy A pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void dgemv_ (char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

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
	
	float dot (int n, float* dx, float* dy, int incx, int incy) {
		return sdot_ (&n, dx, &incx, dy, &incy);
	}
	
	double dot (int n, double* dx, double* dy, int incx, int incy) {
		return ddot_ (&n, dx, &incx, dy, &incy);
	}	
	
	void add_scaled (int n, float da, float *dx, float *dy, int incx, int incy) {
		saxpy_ (&n, &da, dx, &incx, dy, &incy);
	}
	
	void add_scaled (int n, double da, double *dx, double *dy, int incx, int incy) {
		daxpy_ (&n, &da, dx, &incx, dy, &incy);
	}
	
	void matrix_vector_multiply (int m, int n, float alpha, float *a, float *x, float beta, float *y, int lda, int incx, int incy) {
		char charN = 'N';
		
		assert (x != y);
		
		if (lda == -1) {
			lda = m;
		}
				
		sgemv_ (&charN, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
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