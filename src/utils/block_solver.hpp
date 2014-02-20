/*!**********************************************************************
 * \file block_solver.hpp
 * /Users/justinbrown/Dropbox/pisces/src/utils
 * 
 * Created by Justin Brown on 2013-11-17.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BLOCK_SOLVER_HPP_8C3ZNSDI
#define BLOCK_SOLVER_HPP_8C3ZNSDI

namespace utils
{
	void p_block_matrix_factorize (int id, int np, int n, int ntop, int nbot, double* a, int* ipiv, double *x, int *xipiv, int* ns, int *info, int lda = -1, int ldx = -1);
	
	void p_block_matrix_solve (int id, int np, int n, int ntop, int nbot, double* a, int* ipiv, double* b, double *x, int *xipiv, int *ns, int *info, int nrhs = 1, int lda = -1, int ldx = -1, int ldb = -1);
	
	void p_block_tridiag_factorize (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, int* ipiv, double *x, int *xipiv, int *info, int nrhs = 1, int lda = -1);
	
	void p_block_tridiag_solve (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, int* ipiv, double* b, double *x, int *xipiv, int *info, int nrhs = 1, int lda = -1, int ldb = -1);
} /* utils */

#endif /* end of include guard: BLOCK_SOLVER_HPP_8C3ZNSDI */
