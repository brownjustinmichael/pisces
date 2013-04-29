/*!***********************************************************************
 * \file utils.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-29.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef UTILS_HPP_OOQZAQJE
#define UTILS_HPP_OOQZAQJE

#include <cstddef>

namespace utils
{
	void copy (int n, double* dx, double* dy, int incx = 1, int incy = 1);
	
	void scale (int n, double da, double* dx, int incx = 1);
	
	void add_scaled (int n, double da, double *dx, double *dy, int incx = 1, int incy = 1);
	
	void matrix_vector_multiply (int m, int n, double alpha, double *a, double *x, double beta, double *y, int lda = -1, int incx = 1, int incy = 1);
	
	void matrix_factorize (int m, int n, double* a, int *ipiv, int *info = NULL, int lda = -1);
	
	void matrix_solve (int n, double* a, int* ipiv, double* b, int *info = NULL, int nrhs = 1, int lda = -1, int ldb = -1);
} /* utils */

#endif /* end of include guard: UTILS_HPP_OOQZAQJE */
