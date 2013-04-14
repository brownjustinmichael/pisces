/*!***********************************************************************
 * \file solver.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-13.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_HPP_N0BUAX6H
#define SOLVER_HPP_N0BUAX6H

#include <vector>
#include "../plan.hpp"

/*!*******************************************************************
 * \brief Function from BLAS that copies a double array to another in place
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param dx The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 * \param dy The array to which the data are copied
 * \param incy A pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void dcopy_(int *n, double *x, int *incx, double *y, int *incy);

/*!*******************************************************************
 * \brief Function from BLAS that swaps a double array with another in place
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param x The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 * \param y The array to which the data are copied
 * \param incy A pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void dswap_(int *n, double *dx, int *incx, double *dy, int *incy);

extern "C" void daxpy_ (int *n, double *da, double *dx, int *incx, double *y, int *incy);

extern "C" void dscal_ (int *n, double *da, double *dx, int *incx);

/*!*******************************************************************
 * \brief Function from LAPACK that factorizes the matrix a by LU decomposition
 * 
 * \param m A pointer to the number of rows in a
 * \param n A pointer to the number of columns in a
 * \param a A double matrix to be overwritten with its LU decomposition
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

namespace solver
{
	enum flags {
		factorized = 0x04
	};
	
	class solver
	{
	public:
		virtual ~solver () {}
		virtual void solve (int *flags = 0x00) = 0;
	};

	class lapack_solver : public solver
	{
	public:
		lapack_solver (int i_n, double *i_data_in, double *i_rhs, double *i_matrix, double *i_data_out) {
			n = i_n;
			data_in = i_data_in;
			rhs = i_rhs;
			matrix = i_matrix;
			data_out = i_data_out;
			ipiv.resize (i_n, 0);
		};
		virtual ~lapack_solver () {}
		void solve (int *flags = 0x00);

	private:
		int n;
		double *data_in;
		double *rhs;
		double *matrix;
		double *data_out;
		std::vector<int> ipiv;
	};
} /* solver */

#endif /* end of include guard: SOLVER_HPP_N0BUAX6H */
