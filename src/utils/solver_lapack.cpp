/*!**********************************************************************
 * \file solver_lapack.cpp
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "solver_utils.hpp"
#include "lapack.hpp"

#define LAPACK

namespace utils
{
	void matrix_factorize (int m, int n, double* a, int *ipiv, int *info, int lda) {
		if (!info) {
			int iinfo;
			info = &iinfo;
		}
		
		if (lda == -1) {
			lda = m;
		}
		
		dgetrf_ (&m, &n, a, &lda, ipiv, info);
	}
	
	void matrix_solve (int n, double* a, int* ipiv, double* b, int *info, int nrhs, int lda, int ldb) {
		char charN = 'N';
		
		if (!info) {
			int iinfo;
			info = &iinfo;
		}
		
		if (lda == -1) {
			lda = n;
		}
		
		if (ldb == -1) {
			ldb = n;
		}
		
		dgetrs_ (&charN, &n, &nrhs, a, &lda, ipiv, b, &ldb, info);
	}
} /* utils */