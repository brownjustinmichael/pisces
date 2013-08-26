/*!**********************************************************************
 * \file solver_utils_cuda.hpp
 * /Users/justinbrown/Dropbox/spectral_element/src
 * 
 * Created by Justin Brown on 2013-08-22.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_UTILS_HPP_T5OEKE5A
#define SOLVER_UTILS_HPP_T5OEKE5A

namespace cuda
{
	namespace utils
	{
		void matrix_solve (int n, double* a, int* ipiv, double* b, int nrhs = 1, int lda = -1, int ldb = -1);

		void matrix_solve (int n, float* a, int* ipiv, float* b, int nrhs = 1, int lda = -1, int ldb = -1);
	} /* utils */
} /* cuda */

#endif /* end of include guard: SOLVER_UTILS_HPP_T5OEKE5A */
