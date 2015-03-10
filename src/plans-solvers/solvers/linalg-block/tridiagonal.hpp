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
		void tridiag_factorize (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, int* ipiv, double *x, int *xipiv, int *info, int nrhs = 1, int lda = -1);
	
		void tridiag_solve (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, int* ipiv, double* b, double *x, int *xipiv, int *info, int nrhs = 1, int lda = -1, int ldb = -1);
	} /* block */
} /* linalg */

#endif /* end of include guard: BLOCK_TRIDIAGONAL_HPP_BEF6CB65 */
