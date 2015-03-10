/*!***********************************************************************
 * \file utils.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-29.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef UTILS_HPP_OOQZAQJE
#define UTILS_HPP_OOQZAQJE

/*!**********************************************************************
 * \namespace linalg
 * 
 * \brief A namespace containing a set of useful linear algebra functions from BLAS and LAPACK.
 ************************************************************************/
namespace linalg
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
	void copy (int n, const float* x, float* y, int incx = 1, int incy = 1);

	/*!*******************************************************************
	 * \brief Copies a double array into another
	 * 
	 * \param n An integer number of elements to copy
	 * \param x A double array of values to be copied
	 * \param y A double array to be copied to
	 * \param incx The integer spacing of elements in dx
	 * \param incy The integer spacing of elements in dy
	 *********************************************************************/
	void copy (int n, const double* x, double* y, int incx = 1, int incy = 1);
	
	/*!**********************************************************************
	 * \brief Copies a float matrix into another
	 * 
	 * \param n An integer number of elements on the main axis
	 * \param m An integer number of elements on the secondary axis
	 * \param x A float array of values to be copied
	 * \param y A float array to be copied to
	 * \param ldx The leading dimension of x
	 * \param ldy The leading dimension of y
	 ************************************************************************/
	void matrix_copy (int n, int m, float *x, float *y, int ldx = -1, int ldy = -1);
	
	/*!**********************************************************************
	 * \brief Copies a double matrix into another
	 * 
	 * \param n An integer number of elements on the main axis
	 * \param m An integer number of elements on the secondary axis
	 * \param x A double array of values to be copied
	 * \param y A double array to be copied to
	 * \param ldx The leading dimension of x
	 * \param ldy The leading dimension of y
	 ************************************************************************/
	void matrix_copy (int n, int m, double *x, double *y, int ldx = -1, int ldy = -1);

	/*!**********************************************************************
	 * \brief Copies an int matrix into another
	 * 
	 * \param n An integer number of elements on the main axis
	 * \param m An integer number of elements on the secondary axis
	 * \param x A int array of values to be copied
	 * \param y A int array to be copied to
	 * \param ldx The leading dimension of x
	 * \param ldy The leading dimension of y
	 ************************************************************************/
	void matrix_copy (int n, int m, int *x, int *y, int ldx = -1, int ldy = -1);
	
	/*!**********************************************************************
	 * \brief Swap two float matrices
	 * 
	 * \param n An integer number of elements on the main axis
	 * \param m An integer number of elements on the secondary axis
	 * \param x A int array of values to be copied
	 * \param y A int array to be copied to
	 * \param ldx The leading dimension of x
	 * \param ldy The leading dimension of y
	 ************************************************************************/
	void matrix_switch (int n, int m, float *x, float *y, int ldx = -1, int ldy = -1);
	
	/*!**********************************************************************
	 * \brief Swap two double matrices
	 * 
	 * \param n An integer number of elements on the main axis
	 * \param m An integer number of elements on the secondary axis
	 * \param x A int array of values to be copied
	 * \param y A int array to be copied to
	 * \param ldx The leading dimension of x
	 * \param ldy The leading dimension of y
	 ************************************************************************/
	void matrix_switch (int n, int m, double *x, double *y, int ldx = -1, int ldy = -1);

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
	
	/*!**********************************************************************
	 * \brief Scale a float matrix
	 * 
	 * \param n An integer number of elements on the main axis
	 * \param m An integer number of elements on the secondary axis
	 * \param a A float value to scale by
	 * \param x A float array to be scaled
	 * \param ldx The leading dimension of x
	 * \param ldy The leading dimension of y
	 ************************************************************************/
	void matrix_scale (int n, int m, float a, float *x, int ldx = -1);
	
	/*!**********************************************************************
	 * \brief Scale a double matrix
	 * 
	 * \param n An integer number of elements on the main axis
	 * \param m An integer number of elements on the secondary axis
	 * \param a A double value to scale by
	 * \param x A double array to be scaled
	 * \param ldx The leading dimension of x
	 * \param ldy The leading dimension of y
	 ************************************************************************/
	void matrix_scale (int n, int m, double a, double *x, int ldx = -1);
	
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
	
	/*!**********************************************************************
	 * \brief Add a scaled matrix to another float matrix
	 * 
	 * \param n An integer number of elements on the main axis
	 * \param m An integer number of elements on the secondary axis
	 * \param da A float value to scale dx
	 * \param dx A float array of values to be scaled and copied
	 * \param dy A float array to be added to
	 * \param ldx The leading dimension of x
	 * \param ldy The leading dimension of y
	 ************************************************************************/
	void matrix_add_scaled (int n, int m, float da, float *dx, float *dy, int ldx = -1, int ldy = -1);
	
	/*!**********************************************************************
	 * \brief Add a scaled matrix to another double matrix
	 * 
	 * \param n An integer number of elements on the main axis
	 * \param m An integer number of elements on the secondary axis
	 * \param da A double value to scale dx
	 * \param dx A double array of values to be scaled and copied
	 * \param dy A double array to be added to
	 * \param ldx The leading dimension of x
	 * \param ldy The leading dimension of y
	 ************************************************************************/
	void matrix_add_scaled (int n, int m, double da, double *dx, double *dy, int ldx = -1, int ldy = -1);
	
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
	 * \brief Perform the matrix-matrix multiplication y = alpha * a * b + beta * c
	 * 
	 * \param m An integer number of rows in c
	 * \param n An integer number of columns in c
	 * \param k An integer number of columns in a/rows in b
	 * \param alpha The float alpha
	 * \param a The float matrix a (Fortran format)
	 * \param b The float matrix b
	 * \param beta The float beta
	 * \param c The float matrix c
	 * \param lda The integer leading dimension of a
	 * \param ldb The integer leading dimension of b
	 * \param ldc The integer leading dimension of c
	 *********************************************************************/
	void matrix_matrix_multiply (int m, int n, int k, float alpha, float *a, float *b, float beta, float *c, int lda = -1, int ldb = -1, int ldc = -1);

	/*!*******************************************************************
	 * \brief Perform the matrix-matrix multiplication y = alpha * a * b + beta * c
	 * 
	 * \param m An integer number of rows in c
	 * \param n An integer number of columns in c
	 * \param k An integer number of columns in a/rows in b
	 * \param alpha The double alpha
	 * \param a The double matrix a (Fortran format)
	 * \param b The double matrix b
	 * \param beta The double beta
	 * \param c The double matrix c
	 * \param lda The integer leading dimension of a
	 * \param ldb The integer leading dimension of b
	 * \param ldc The integer leading dimension of c
	 *********************************************************************/
	void matrix_matrix_multiply (int m, int n, int k, double alpha, double *a, double *b, double beta, double *c, int lda = -1, int ldb = -1, int ldc = -1);
	
	/*!*******************************************************************
	 * \brief Perform the diagonal matrix-vector multiplication y = alpha * a * x + beta * y
	 * 
	 * \param m An integer number of rows in a
	 * \param n An integer number of columns in a
	 * \param alpha The float alpha
	 * \param a The float matrix a (Fortran format)
	 * \param x The float array x
	 * \param beta The float beta
	 * \param y The float array y
	 * \param inca The integer spacing of elements in a
	 * \param incx The integer spacing of elements in x
	 * \param incy The integer spacing of elements in y
	 *********************************************************************/
	void diagonal_multiply (int n, float alpha, float *a, float *x, float beta, float *y, int inca = 1, int incx = 1, int incy = 1);
	
	/*!*******************************************************************
	 * \brief Perform the diagonal matrix-vector multiplication y = alpha * a * x + beta * y
	 * 
	 * \param m An integer number of rows in a
	 * \param n An integer number of columns in a
	 * \param alpha The double alpha
	 * \param a The double matrix a (Fortran format)
	 * \param x The double array x
	 * \param beta The double beta
	 * \param y The double array y
	 * \param inca The integer spacing of elements in a
	 * \param incx The integer spacing of elements in x
	 * \param incy The integer spacing of elements in y
	 *********************************************************************/
	void diagonal_multiply (int n, double alpha, double *a, double *x, double beta, double *y, int inca = 1, int incx = 1, int incy = 1);

} /* linalg */

#endif /* end of include guard: UTILS_HPP_OOQZAQJE */
