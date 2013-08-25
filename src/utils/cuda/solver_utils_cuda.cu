/*!**********************************************************************
 * \file solver_utils_cuda.cu
 * /Users/justinbrown/Dropbox/spectral_element/src
 * 
 * Created by Justin Brown on 2013-08-22.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "solver_utils_cuda.hpp"

namespace cuda
{
	namespace utils
	{
		void matrix_solve (int n, double* a, int* ipiv, double* b, int nrhs = 1, int lda = -1, int ldb = -1) {
		    if (n == 0 || nrhs == 0) {
			return;
		    }
			
			if (lda == -1) {
				lda = n;
			}
			if (ldb == -1) {
				ldb = n;
			}
			
		    if (n < 0) {
				throw -2;
		    } else if (nrhs < 0) {
				throw -3;
		    } else if (lda < max(1,n)) {
				throw -5;
		    } else if (ldb < max(1,n)) {
				throw -8;
		    }
	
			swap <<<1, min (nrhs, 256)>>> (nrhs, b, ldb, 1, n, ipiv, 1);
			operation_left_lower <<<min (nrhs, 1024), min (n, 256)>>> (n, nrhs, a, lda, b, ldb, false);
			operation_left_upper <<<min (nrhs, 1024), min (n, 256)>>> (n, nrhs, a, lda, b, ldb, true);
		}
	} /* utils */
} /* cuda */