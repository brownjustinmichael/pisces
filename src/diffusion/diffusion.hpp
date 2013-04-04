// 
//! \file diffusion.hpp
//  diffusion
// 
//  This is the header file for the basic diffusion classes.
//  
//  Created by Justin Brown on 2013-03-22.
//  Copyright 2013 Justin Brown. All rights reserved.
// 

#ifndef diffusion_H_1AIIK1RA
#define diffusion_H_1AIIK1RA

#include <vector>
#include "../plan.hpp"
#include "../collocation/collocation.hpp"

//! \brief Function from blas that copies a double array to another in place
//
//! \param n a pointer to an integer number of elements in x to copy to y
//! \param x the array from which the data are copied
//! \param incx a pointer to an integer spacing of elements in x
//! \param y the array to which the data are copied
//! \param incy a pointer to an integer spacing of elements in y
extern "C" void   dcopy_(int *n, double *x, int *incx, double *y, int *incy);

//! \brief Function from blas that solves a packed triangular matrix equation
//
//! \param uplo a pointer to a character, "U" for upper triangular matrix, "L" for lower triangular matrix
//! \param trans (no idea, just set this to a pointer to "N")
//! \param diag a pointer to a character, "U" for unit diagonal, "N" for anything else
//! \param n a pointer to an integer number of elements in x to copy to y
//! \param a the triangular packed array of matrix elements
//! \param x the array of right-hand-side elements in the equation which will be replaced with the resulting array
//! \param incx a pointer to an integer spacing of elements in x
extern "C" void   dtpsv_(char *uplo, char *trans, char *diag, int *n, double *a, double *x, int *incx);

extern "C" void dgtsv_(int *n, int *nrhs, double *dl, double *d, double *du, double *b, int *ldb, int *info);

extern "C" void dgesv_ (int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

extern "C" void dgemv_ (char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

extern "C" void dgbsv_(int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info);

namespace diffusion
{
	//! \brief Subclass of operation, one implementation of diffusion for data expressed as a sum of Chebyshev polynomials in 1D
	//
	//! This implementation is a simple forward time-difference scheme, solving only for the diffusion component of the equations. The class solves this implicitly with two matrix equations, one for the even and one for the odd Chebyshev polynomials. It uses the BLAS library to do so.
	class cheb_1D : public plan
	{
	private:
		double coeff; //!< a double that represents the coefficient in front of the diffusion term in the differential equation
		double alpha;
		int n; //!< an integer number of data elements (grid points) that cheb_1D will be built to handle
		double *data_in; //!< a double pointer to the input data
		double *data_out; //!< a double pointer to the output data; if data_in == data_out, the operation is done in place
		std::vector<double> diffusion_matrix; //!< a 1D vector to be filled with the matrix equation for the Chebyshev polynomials
		std::vector<double> pre_matrix;
		std::vector<double> temp;
		std::vector<int> ipiv;
		collocation::cheb_grid *cheb;
	public:
		//! \param i_coeff a double containing the coefficient in front of the diffusion term in the differential equation
		//! \param i_n an integer number of data elements (grid points) that cheb_1D will be built to tackle
		//! \param i_data_in a double pointer pointing to the input data
		//! \param i_data_out a double pointer pointing to the output data
		cheb_1D (double i_coeff, double alpha, int i_n, double *i_data_in, double *i_data_out = NULL);
		virtual ~cheb_1D () {delete cheb;}
		//! \brief Execute the operation on the data for a given timestep duration
		//
		//! \param timestep a double duration over which the diffusion step will happen
		void execute (double timestep);
	};
	
	class fixed_angle_1D : public plan
	{
	private:
		double coeff;
		int n;
		double *data_in;
		double *data_out;
		std::vector<double> diagonal_left;
		std::vector<double> diagonal_right;
		std::vector<double> diagonal;
	public:
		fixed_angle_1D (double i_coeff, int i_n, double *i_data_in, double *i_data_out = NULL);
		void execute (double timestep);
	};
} /* diffusion */

#endif /* end of include guard: diffusion_H_1AIIK1RA */
