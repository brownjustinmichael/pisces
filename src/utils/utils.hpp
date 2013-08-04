/*!***********************************************************************
 * \file utils.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-29.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef UTILS_HPP_OOQZAQJE
#define UTILS_HPP_OOQZAQJE

#include <cstddef>

namespace utils
{
	/*!*******************************************************************
	 * \brief Copies a double array into another
	 * 
	 * \param n An integer number of elements to copy
	 * \param dx A double array of values to be copied
	 * \param dy A double array to be copied to
	 * \param incx The integer spacing of elements in dx
	 * \param incy The integer spacing of elements in dy
	 *********************************************************************/
	void copy (int n, double* dx, double* dy, int incx = 1, int incy = 1);
	
	/*!*******************************************************************
	 * \brief Scales a double array by a constant
	 * 
	 * \param n An integer number of elements to scale
	 * \param da A double to scale the array by
	 * \param dx A double array to scale
	 * \param incx The integer spacing of elements in dx
	 *********************************************************************/
	void scale (int n, double da, double* dx, int incx = 1);

	/*!*******************************************************************
	 * \brief Takes the dot product of two arrays
	 * 
	 * \param n An integer number of elements to dot
	 * \param dx The first double array
	 * \param dy The second double array
	 * \param incx The integer spacing of elements in dx
	 * \param incy The integer spacing of elements in dy
	 * 
	 * \return The double dot product
	 *********************************************************************/
	double dot (int n, double* dx, double* dy, int incx = 1, int incy = 1);
	
	/*!*******************************************************************
	 * \brief Perform the operation dy = da * dx + dy
	 * 
	 * \param n An integer number of elements to operate on
	 * \param da The double da
	 * \param dx The double array dx
	 * \param dy The double array dy which is replaced with the output
	 * \param incx The integer spacing of elements in dx
	 * \param incy The integer spacing of elements in dy
	 *********************************************************************/
	void add_scaled (int n, double da, double *dx, double *dy, int incx = 1, int incy = 1);
	
	/*!*******************************************************************
	 * \brief Perform the matrix-vector multiplication y = alpha * a * x + beta * y
	 * 
	 * \param m An integer number of rows in a
	 * \param n An integer number of columns in a
	 * \param alpha The double alpha
	 * \param a The double matrix a (Fortran format)
	 * \param x The double array x
	 * \param beta The double beta
	 * \param y The double array y
	 * \param lda The integer leading dimension of a
	 * \param incx The integer spacing of elements in x
	 * \param incy The integer spacing of elements in y
	 *********************************************************************/
	void matrix_vector_multiply (int m, int n, double alpha, double *a, double *x, double beta, double *y, int lda = -1, int incx = 1, int incy = 1);
	
	/*!*******************************************************************
	 * \brief LU Decompose the matrix a
	 * 
	 * \param m An integer number of rows in a
	 * \param n An integer number of columns in a
	 * \param a The double matrix a (Fortran format)
	 * \param ipiv An integer array to contain output pivot indices
	 * \param info An integer pointer that points to 0 on success, something else on failure
	 * \param lda The leading dimension of a
	 *********************************************************************/
	void matrix_factorize (int m, int n, double* a, int *ipiv, int *info = NULL, int lda = -1);
	
	/*!*******************************************************************
	 * \brief Solve b = A * x and output result in b
	 * 
	 * \param n An integer number of elements in b
	 * \param a The double matrix a (LU Decomposed)
	 * \param ipiv The integer array of pivot indices from the LU Decomposition
	 * \param b The double array b
	 * \param info An integer pointer that points to 0 on success, something else on failure
	 * \param nrhs The number of right hand sides in b
	 * \param lda The leading dimension of a
	 * \param ldb The leading dimension of b
	 *********************************************************************/
	void matrix_solve (int n, double* a, int* ipiv, double* b, int *info = NULL, int nrhs = 1, int lda = -1, int ldb = -1);

	double interpolate (int n, double* dx, double* dy, double x);
	
	void matrix_interpolate (int n, double* dx, int m, double* dy, int incy, double* douty, double x);
} /* utils */

#endif /* end of include guard: UTILS_HPP_OOQZAQJE */
