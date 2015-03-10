/*!**********************************************************************
 * \file linalg.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef UTILS_SOLVER_HPP_YSBJBB1J
#define UTILS_SOLVER_HPP_YSBJBB1J

#include <cstddef>
#include "mpi/messenger.hpp"

namespace linalg
{
	/*!**********************************************************************
	 * \brief Solve the equation a * b = b, where a is a diagonal matrix
	 * 
	 * \param n The number of rows/columns in a
	 * \param a A pointer to the float matrix a
	 * \param b A pointer to the float vector b
	 * \param inca The integer spacing of a
	 * \param incb The integer spacing of b
	 ************************************************************************/
	void diagonal_solve (int n, float *a, float *b, int inca = 1, int incb = 1);
	
	/*!**********************************************************************
	 * \brief Solve the equation a * b = b, where a is a diagonal matrix
	 * 
	 * \param n The number of rows/columns in a
	 * \param a A pointer to the double matrix a
	 * \param b A pointer to the double vector b
	 * \param inca The integer spacing of a
	 * \param incb The integer spacing of b
	 ************************************************************************/
	void diagonal_solve (int n, double *a, double *b, int inca = 1, int incb = 1);
	
	/*!*******************************************************************
	 * \brief LU Decompose the matrix a
	 * 
	 * \param m An integer number of rows in a
	 * \param n An integer number of columns in a
	 * \param a The float matrix a (Fortran format)
	 * \param ipiv An integer array to contain output pivot indices
	 * \param info An integer pointer that points to 0 on success, something else on failure
	 * \param lda The leading dimension of a
	 *********************************************************************/
	void matrix_factorize (int m, int n, float* a, int *ipiv, int *info = NULL, int lda = -1);
	
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
	 * \param a The float matrix a (LU Decomposed)
	 * \param ipiv The integer array of pivot indices from the LU Decomposition
	 * \param b The float array b
	 * \param info An integer pointer that points to 0 on success, something else on failure
	 * \param nrhs The number of right hand sides in b
	 * \param lda The leading dimension of a
	 * \param ldb The leading dimension of b
	 *********************************************************************/
	void matrix_solve (int n, float* a, int* ipiv, float* b, int *info = NULL, int nrhs = 1, int lda = -1, int ldb = -1);

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
	
	/*!**********************************************************************
	 * \brief Factorize a banded matrix
	 * 
	 * \param m The integer number of rows in a
	 * \param n The integer number of columns in a
	 * \param kl The integer number of rows below the diagonal
	 * \param ku The integer number of rows above the diagonal
	 * \param a The double matrix to factorize
	 * \param ipiv The integer array of pivot indices from the decomposition
	 * \param info An integer pointer that points to 0 on success, something else on failure
	 * \param lda The integer leading dimension of a (must be greater than 2 * kl + ku + 1)
	 * 
	 * Warning: the matrix a must be an m * (2 * kl + ku + 1) matrix to fit the decomposition
	 ************************************************************************/
	void matrix_banded_factorize (int m, int n, int kl, int ku, double* a, int *ipiv, int *info = NULL, int lda = -1);
	
	/*!**********************************************************************
	 * \brief Solve a factorized banded matrix
	 * 
	 * \param n The integer number of rows/columns in a
	 * \param kl The integer number of rows below the diagonal
	 * \param ku The integer number of rows above the diagonal
	 * \param a The double factorized matrix
	 * \param ipiv The integer array of pivot indices from the decomposition
	 * \param b The double right hand side array/output on result
	 * \param info An integer pointer that points to 0 on success, something else on failure
	 * \param nrhs The integer number of right hand sides in b
	 * \param lda The integer leading dimension of a (must be greater than 2 * kl + ku + 1)
	 * \param ldb The integer leading dimension of b
	 * 
	 * Warning: the matrix a must be an m * (2 * kl + ku + 1) matrix to fit the decomposition
	 ************************************************************************/
	void matrix_banded_solve (int n, int kl, int ku, double* a, int* ipiv, double* b, int *info = NULL, int nrhs = 1, int lda = -1, int ldb = -1);
	
	/*!**********************************************************************
	 * \brief Factorize a tridiagonal matrix
	 * 
	 * \param n The integer number of columns in the matrix
	 * \param dl The float vector of subdiagonal elements
	 * \param d The float vector of diagonal elements
	 * \param du The float vector of superdiagonal elements
	 * \param du2 The float vector of supersuperdiagonal elements
	 * \param ipiv The integer array of pivot indices from the decomposition
	 * \param info An integer pointer that points to 0 on success, something else on failure
	 ************************************************************************/
	void tridiagonal_factorize (int n, float *dl, float *d, float *du, float *du2, int *ipiv, int *info = NULL);
	
	/*!**********************************************************************
	 * \brief Factorize a tridiagonal matrix
	 * 
	 * \param n The integer number of columns in the matrix
	 * \param dl The double vector of subdiagonal elements
	 * \param d The double vector of diagonal elements
	 * \param du The double vector of superdiagonal elements
	 * \param du2 The double vector of supersuperdiagonal elements
	 * \param ipiv The integer array of pivot indices from the decomposition
	 * \param info An integer pointer that points to 0 on success, something else on failure
	 ************************************************************************/
	void tridiagonal_factorize (int n, double *dl, double *d, double *du, double *du2, int *ipiv, int *info = NULL);
	
	/*!**********************************************************************
	 * \brief Solve a tridiagonal matrix
	 * 
	 * \param n The integer number of columns in the matrix
	 * \param dl The float vector of subdiagonal elements
	 * \param d The float vector of diagonal elements
	 * \param du The float vector of superdiagonal elements
	 * \param du2 The float vector of supersuperdiagonal elements
	 * \param ipiv The integer array of pivot indices from the decomposition
	 * \param b The float array right hand side to be overwritten
	 * \param info An integer pointer that points to 0 on success, something else on failure
	 * \param nrhs The integer number of right hand sides
	 * \param ldb The integer leading dimension of b
	 ************************************************************************/
	void tridiagonal_solve (int n, float *dl, float *d, float *du, float *du2, int *ipiv, float *b, int *info = NULL, int nrhs = 1, int ldb = -1);
	
	/*!**********************************************************************
	 * \brief Solve a tridiagonal matrix
	 * 
	 * \param n The integer number of columns in the matrix
	 * \param dl The double vector of subdiagonal elements
	 * \param d The double vector of diagonal elements
	 * \param du The double vector of superdiagonal elements
	 * \param du2 The double vector of supersuperdiagonal elements
	 * \param ipiv The integer array of pivot indices from the decomposition
	 * \param b The double array right hand side to be overwritten
	 * \param info An integer pointer that points to 0 on success, something else on failure
	 * \param nrhs The integer number of right hand sides
	 * \param ldb The integer leading dimension of b
	 ************************************************************************/
	void tridiagonal_solve (int n, double *dl, double *d, double *du, double *du2, int *ipiv, double *b, int *info = NULL, int nrhs = 1, int ldb = -1);
	
	/*!**********************************************************************
	 * \brief Solve a tridiagonal matrix directly
	 * 
	 * \param n The integer number of rows/columns in the matrix
	 * \param sub The float vector of subdiagonal elements
	 * \param diag The float vector of diagonal elements
	 * \param sup The float vector of superdiagonal elements
	 * \param b The float right hand side array/output on result
	 * \param info An integer pointer that points to 0 on success, something else on failure
	 * \param nrhs The integer number of right hand sides in b
	 * \param lda The integer leading dimension of a
	 * \param ldb The integer leading dimension of b
	 * 
	 * Warning: this function overwrites the original matrices and can't be used again until the matrices are reset
	 ************************************************************************/
	void tridiagonal_solve (int n, float *sub, float *diag, float *sup, float *b, int *info = NULL, int nrhs = 1, int ldb = -1);
	
	/*!**********************************************************************
	 * \brief Solve a tridiagonal matrix directly
	 * 
	 * \param n The integer number of rows/columns in the matrix
	 * \param sub The double vector of subdiagonal elements
	 * \param diag The double vector of diagonal elements
	 * \param sup The double vector of superdiagonal elements
	 * \param b The double right hand side array/output on result
	 * \param info An integer pointer that points to 0 on success, something else on failure
	 * \param nrhs The integer number of right hand sides in b
	 * \param lda The integer leading dimension of a
	 * \param ldb The integer leading dimension of b
	 * 
	 * Warning: this function overwrites the original matrices and can't be used again until the matrices are reset
	 ************************************************************************/
	void tridiagonal_solve (int n, double *sub, double *diag, double *sup, double *b, int *info = NULL, int nrhs = 1, int ldb = -1);
} /* linalg */

#endif /* end of include guard: UTILS_SOLVER_HPP_YSBJBB1J */
