/*!**********************************************************************
 * \file block_solver.hpp
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
		void banded_factorize (int id, int np, int n, int kl, int ku, double* matrix, int* ipiv, double *x, int *xipiv, double *bufferl, double *bufferr, int *info = NULL, int nrhs = 1, int lda = -1, int ldaa = -1);
	
		void banded_solve (int id, int np, int n, int kl, int ku, double* matrix, int* ipiv, double* b, double *x, int *xipiv, double *bufferl, double *bufferr, int *info = NULL, int nrhs = 1, int lda = -1, int ldaa = -1, int ldb = -1);
	} /* block */
} /* linalg */

#endif /* end of include guard: BLOCK_SOLVER_HPP_8C3ZNSDI */
