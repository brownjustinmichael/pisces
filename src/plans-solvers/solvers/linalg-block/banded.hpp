/*!**********************************************************************
 * \file banded.hpp
 * /Users/justinbrown/Dropbox/pisces/src/utils
 * 
 * Created by Justin Brown on 2013-11-17.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BLOCK_SOLVER_HPP_8C3ZNSDI
#define BLOCK_SOLVER_HPP_8C3ZNSDI

namespace linalg
{
	namespace block
	{
		/*!**********************************************************************
		 * \brief Factorize a banded block matrix
		 * 
		 * \param id The mpi id of this process
		 * \param np The number of mpi processes/number of blocks in the matrix
		 * \param n The number of elements in the non-overlapping part of the matrix
		 * \param kl The number of subdiagonal elements
		 * \param ku The number of superdiagonal elements
		 * \param matrix The matrix that holds the data (in banded LAPACK format, must be n * (2 * kl + ku + 1 + ntop + nbot))
		 * \param ipiv The integer array to hold the LAPACK swap information
		 * \param x Some additional memory needed for the solve, must be (ntop + nbot)^2 on all but process 0, for which it must be (sum (ntop + nbot))^2 over all processes
		 * \param xipiv The integer array to hold the LAPACK swap information for x
		 * \param bufferl A buffer that is nrhs * kl * n
		 * \param bufferr A buffer that is nrhs * ku * n
		 * \param info A pointer to the result state of the solve, 0 if successful
		 * \param nrhs The number of right hand sides in the matrix
		 * \param lda The leading dimension of matrix, if < 0, 2 * kl + ku + 1
		 * \param ldaa The next leading dimension of matrix, if < 0, n + ku + kl + ntop + nbot
		 * 
		 * The input matrix should be in banded LAPACK format. To handle the overlap region, which wouldn't normally appear in such a format, extend the matrix by kl on the left and ku on the right. The overlapping data should be included in the same format. For details on the algorithm, see linalg::block::matrix_factorize.
		 ************************************************************************/
		void banded_factorize (int id, int np, int n, int kl, int ku, double* matrix, int* ipiv, double *x, int *xipiv, double *bufferl, double *bufferr, double *buffer, int *info = NULL, int nrhs = 1, int lda = -1, int ldaa = -1);
		
		/*!**********************************************************************
		 * \brief Solve a banded block matrix
		 * 
		 * \param id The mpi id of this process
		 * \param np The number of mpi processes/number of blocks in the matrix
		 * \param n The number of elements in the non-overlapping part of the matrix
		 * \param kl The number of subdiagonal elements
		 * \param ku The number of superdiagonal elements
		 * \param matrix The matrix that holds the factorized data (in banded LAPACK format, must be n * (2 * kl + ku + 1 + ntop + nbot))
		 * \param ipiv The integer array to hold the LAPACK swap information
		 * \pararm b The right hand side data
		 * \param x Some additional memory needed for the solve, must be (ntop + nbot)^2 on all but process 0, for which it must be (sum (ntop + nbot))^2 over all processes
		 * \param xipiv The integer array to hold the LAPACK swap information for x
		 * \param bufferl A buffer that is nrhs * kl * n from factorize
		 * \param bufferr A buffer that is nrhs * ku * n from factorize
		 * \param info A pointer to the result state of the solve, 0 if successful
		 * \param nrhs The number of right hand sides in the matrix
		 * \param lda The leading dimension of matrix, if < 0, 2 * kl + ku + 1
		 * \param ldaa The next leading dimension of matrix, if < 0, n + ku + kl + ntop + nbot
		 * \param ldb The leading dimension of b, if < 0, n + ntop + nbot
		 * 
		 * For details on the algorithm, see linalg::block::matrix_solve
		 ************************************************************************/
		void banded_solve (int id, int np, int n, int kl, int ku, double* matrix, int* ipiv, double* b, double *x, int *xipiv, double *bufferl, double *bufferr, int *info = NULL, int nrhs = 1, int lda = -1, int ldaa = -1, int ldb = -1);
	} /* block */
} /* linalg */

#endif /* end of include guard: BLOCK_SOLVER_HPP_8C3ZNSDI */
