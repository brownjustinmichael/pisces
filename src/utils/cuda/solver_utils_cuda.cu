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
		void matrix_solve (int n, double* a, int* ipiv, double* b, int *info = NULL, int nrhs = 1, int lda = -1, int ldb = -1) {
			
		}
	} /* utils */
} /* cuda */