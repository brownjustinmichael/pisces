/*!***********************************************************************
 * \file lapack.hpp
 * Spectral Element
 * 
 * This file simply contains the LAPACK and BLAS routines that can be used throughout the code.
 * 
 * Created by Justin Brown on 2013-04-29.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef LAPACK_HPP_MEIP1OZM
#define LAPACK_HPP_MEIP1OZM

/*!*******************************************************************
 * \brief Function from BLAS that calculates a dot product
 * 
 * \param n A pointer to an integer number of elements in dx and dy to dot
 * \param dx The first double array of data
 * \param incx A pointer to an integer spacing of elements in dx
 * \param dy The second double array of data
 * \param incy A pointer to an integer spacing of elements in dy
 *********************************************************************/
extern "C" double ddot_ (int *n, double *x, int *incx, double *y, int *incy);

/*!*******************************************************************
 * \brief Function from BLAS that copies a float array to another in place
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param dx The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 * \param dy The array to which the data are copied
 * \param incy A pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void scopy_ (int *n, float *x, int *incx, float *y, int *incy);

/*!*******************************************************************
 * \brief Function from BLAS that copies a double array to another in place
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param dx The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 * \param dy The array to which the data are copied
 * \param incy A pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void dcopy_ (int *n, double *x, int *incx, double *y, int *incy);

/*!*******************************************************************
 * \brief Function from BLAS that scales a float array
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param da The float by which to scale the data
 * \param dx The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 *********************************************************************/
extern "C" void sscal_ (int *n, float *a, float *x, int *incx);

/*!*******************************************************************
 * \brief Function from BLAS that scales a double array
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param da The double by which to scale the data
 * \param dx The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 *********************************************************************/
extern "C" void dscal_ (int *n, double *a, double *x, int *incx);

/*!*******************************************************************
 * \brief Function from BLAS for vector-vector addition (dy = da * dx + dy)
 * 
 * \param n A pointer to the integer number of values in dy
 * \param da A pointer to the double da
 * \param dx The double vector dx
 * \param incx A pointer to the integer spacing of elements in dx
 * \param dy The double vector dy
 * \param incy A pointer to the integer spacing of elements in dy
 *********************************************************************/
extern "C" void daxpy_ (int *n, double *da, double *dx, int *incx, double *dy, int *incy);

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

/*!*******************************************************************
 * \brief Function from LAPACK that factorizes the matrix a by LU decomposition
 * 
 * \param m A pointer to the number of rows in a
 * \param n A pointer to the number of columns in a
 * \param a A double matrix to be overwritten with its LU decomposition
 * \param lda A double pointer to the integer leading dimension of a
 * \param ipiv An integer array to contain the pivot indices
 * \param info A pointer to an integer indicating success (0 for successful exit)
 *********************************************************************/
extern "C" void dgetrf_ (int *m, int *n, double *a, int *lda, int *ipiv, int *info);

/*!*******************************************************************
 * \brief Function from LAPACK that solves a factorized matrix equation
 * 
 * \param trans A pointer to transposition character ("N" for not transposed, "T" for transposed)
 * \param n A pointer to the number of columns in a
 * \param nrhs A pointer to the number of right hand sides
 * \param a A double matrix to be overwritten with its LU decomposition
 * \param lda A pointer to the integer number of leading dimension of a
 * \param ipiv An integer array to contain the pivot indices
 * \param b The double right hand side array, overwritten with solution
 * \param ldb A pointer to the integer number of leading dimension of b
 * \param info A pointer to an integer indicating success (0 for successful exit)
 *********************************************************************/
extern "C" void dgetrs_ (char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

#endif /* end of include guard: LAPACK_HPP_MEIP1OZM */
