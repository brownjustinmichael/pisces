/*!**********************************************************************
 * \file utils_solver.hpp
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef UTILS_SOLVER_HPP_YSBJBB1J
#define UTILS_SOLVER_HPP_YSBJBB1J

#include <cstddef>

namespace utils
{
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
	
} /* utils */

#endif /* end of include guard: UTILS_SOLVER_HPP_YSBJBB1J */
