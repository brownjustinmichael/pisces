/*!***********************************************************************
 * \file utils.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-29.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef UTILS_HPP_OOQZAQJE
#define UTILS_HPP_OOQZAQJE

namespace utils
{
	/*!*******************************************************************
	 * \brief Copies a double array into another
	 * 
	 * \param n An integer number of elements to copy
	 * \param x A float array of values to be copied
	 * \param y A float array to be copied to
	 * \param incx The integer spacing of elements in dx
	 * \param incy The integer spacing of elements in dy
	 *********************************************************************/
	void copy (int n, float* x, float* y, int incx = 1, int incy = 1);

	/*!*******************************************************************
	 * \brief Copies a double array into another
	 * 
	 * \param n An integer number of elements to copy
	 * \param x A double array of values to be copied
	 * \param y A double array to be copied to
	 * \param incx The integer spacing of elements in dx
	 * \param incy The integer spacing of elements in dy
	 *********************************************************************/
	void copy (int n, double* x, double* y, int incx = 1, int incy = 1);
	
	void matrix_copy (int n, int m, float *x, float *y, int ldx, int ldy);
	
	void matrix_copy (int n, int m, double *x, double *y, int ldx, int ldy);
	
	/*!*******************************************************************
	 * \brief Scales a float array by a constant
	 * 
	 * \param n An integer number of elements to scale
	 * \param a A float to scale the array by
	 * \param x A float array to scale
	 * \param incx The integer spacing of elements in dx
	 *********************************************************************/
	void scale (int n, int a, int* x, int incx = 1);

	/*!*******************************************************************
	 * \brief Scales a float array by a constant
	 * 
	 * \param n An integer number of elements to scale
	 * \param a A float to scale the array by
	 * \param x A float array to scale
	 * \param incx The integer spacing of elements in dx
	 *********************************************************************/
	void scale (int n, float a, float* x, int incx = 1);
	
	/*!*******************************************************************
	 * \brief Scales a double array by a constant
	 * 
	 * \param n An integer number of elements to scale
	 * \param a A double to scale the array by
	 * \param x A double array to scale
	 * \param incx The integer spacing of elements in dx
	 *********************************************************************/
	void scale (int n, double a, double* x, int incx = 1);

	/*!*******************************************************************
	 * \brief Takes the dot product of two arrays
	 * 
	 * \param n An integer number of elements to dot
	 * \param x The first double array
	 * \param y The second double array
	 * \param incx The integer spacing of elements in dx
	 * \param incy The integer spacing of elements in dy
	 * 
	 * \return The double dot product
	 *********************************************************************/
	float dot (int n, float* x, float* y, int incx = 1, int incy = 1);
	
	/*!*******************************************************************
	 * \brief Takes the dot product of two arrays
	 * 
	 * \param n An integer number of elements to dot
	 * \param x The first double array
	 * \param y The second double array
	 * \param incx The integer spacing of elements in dx
	 * \param incy The integer spacing of elements in dy
	 * 
	 * \return The double dot product
	 *********************************************************************/
	double dot (int n, double* x, double* y, int incx = 1, int incy = 1);
	
	/*!*******************************************************************
	 * \brief Perform the operation dy = da * dx + dy
	 * 
	 * \param n An integer number of elements to operate on
	 * \param da The double da
	 * \param x The double array dx
	 * \param y The double array dy which is replaced with the output
	 * \param incx The integer spacing of elements in dx
	 * \param incy The integer spacing of elements in dy
	 *********************************************************************/
	void add_scaled (int n, float da, float *x, float *y, int incx = 1, int incy = 1);

	/*!*******************************************************************
	 * \brief Perform the operation dy = da * dx + dy
	 * 
	 * \param n An integer number of elements to operate on
	 * \param da The double da
	 * \param x The double array dx
	 * \param y The double array dy which is replaced with the output
	 * \param incx The integer spacing of elements in dx
	 * \param incy The integer spacing of elements in dy
	 *********************************************************************/
	void add_scaled (int n, double da, double *x, double *y, int incx = 1, int incy = 1);
	
	void matrix_add_scaled (int n, int m, float da, float *dx, float *dy, int ldx, int ldy);
	
	void matrix_add_scaled (int n, int m, double da, double *dx, double *dy, int ldx, int ldy);
	
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
	void matrix_vector_multiply (int m, int n, float alpha, float *a, float *x, float beta, float *y, int lda = -1, int incx = 1, int incy = 1);
	
	void matrix_matrix_multiply (int m, int n, int k, float alpha, float *a, float *b, float beta, float *c, int lda = -1, int ldb = -1, int ldc = -1);

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

	void matrix_matrix_multiply (int m, int n, int k, double alpha, double *a, double *b, double beta, double *c, int lda = -1, int ldb = -1, int ldc = -1);

	void diagonal_multiply (int n, float alpha, float *a, float *x, float beta, float *y, int inca = 1, int incx = 1, int incy = 1);

	void diagonal_multiply (int n, double alpha, double *a, double *x, double beta, double *y, int inca = 1, int incx = 1, int incy = 1);

} /* utils */

#endif /* end of include guard: UTILS_HPP_OOQZAQJE */
