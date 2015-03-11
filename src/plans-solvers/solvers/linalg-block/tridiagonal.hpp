/*!**********************************************************************
 * \file tridiagonal.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BLOCK_TRIDIAGONAL_HPP_BEF6CB65
#define BLOCK_TRIDIAGONAL_HPP_BEF6CB65

namespace linalg
{
	namespace block
	{
		/*!**********************************************************************
		 * \brief Factorize a tridiagonal block matrix
		 * 
		 * \param id The mpi id of this process
		 * \param np The number of mpi processes/number of blocks in the matrix
		 * \param n The number of elements in the non-overlapping part of the matrix
		 * \param sub The subdiagonal elements in the matrix
		 * \param diag The diagonal elements in the matrix
		 * \param sup The superdiagonal elements in the matrix
		 * \param supsup The supersuperdiagonal elements in the matrix
		 * \param ipiv The integer array to hold the LAPACK swap information
		 * \param x Some additional memory needed for the solve, must be (1 + 1)^2 on all but process 0, for which it must be (sum (1 + 1))^2 over all processes
		 * \param info A pointer to the result state of the solve, 0 if successful
		 * \param nrhs The number of right hand sides in the tridiagonal solve
		 * \param lda The leading dimension of a, if < 0, n + ntop + nbot
		 * 
		 * For more info on the process, see the body of linalg::block::matrix_factorize
		 ************************************************************************/
		void tridiag_factorize (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, int* ipiv, double *x, int *xipiv, int *info, int nrhs = 1, int lda = -1);
		
		/*!**********************************************************************
		 * \brief Solve a tridiagonal block matrix
		 * 
		 * \param id The mpi id of this process
		 * \param np The number of mpi processes/number of blocks in the matrix
		 * \param n The number of elements in the non-overlapping part of the matrix
		 * \param sub The subdiagonal elements in the matrix
		 * \param diag The diagonal elements in the matrix
		 * \param sup The superdiagonal elements in the matrix
		 * \param supsup The supersuperdiagonal elements in the matrix
		 * \param ipiv The integer array to hold the LAPACK swap information
		 * \param b The right hand side to solve (overwritten with solution)
		 * \param x Some additional memory needed for the solve, must be (1 + 1)^2 on all but process 0, for which it must be (sum (1 + 1))^2 over all processes
		 * \param info A pointer to the result state of the solve, 0 if successful
		 * \param nrhs The number of right hand sides in the tridiagonal solve
		 * \param lda The leading dimension of a, if < 0, n + ntop + nbot
		 * \param ldb The leading dimension of b, if < 0, n + ntop + nbot
		 ************************************************************************/
		void tridiag_solve (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, int* ipiv, double* b, double *x, int *xipiv, int *info, int nrhs = 1, int lda = -1, int ldb = -1);
	} /* block */
} /* linalg */

#endif /* end of include guard: BLOCK_TRIDIAGONAL_HPP_BEF6CB65 */
