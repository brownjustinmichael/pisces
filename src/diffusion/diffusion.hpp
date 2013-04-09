/*!***********************************************************************
 * \file diffusion.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef diffusion_H_1AIIK1RA
#define diffusion_H_1AIIK1RA

#include <memory>
#include <vector>
#include "../plan.hpp"
#include "../collocation/collocation.hpp"

/*!*******************************************************************
 * \brief Function from BLAS that copies a double array to another in place
 * 
 * \param n a pointer to an integer number of elements in x to copy to y
 * \param x the array from which the data are copied
 * \param incx a pointer to an integer spacing of elements in x
 * \param y the array to which the data are copied
 * \param incy a pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void   dcopy_(int *n, double *x, int *incx, double *y, int *incy);

extern "C" void dgesv_ (int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

extern "C" void dgemv_ (char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

namespace diffusion
{
	/*!*******************************************************************
	 * \brief Subclass of operation, one implementation of diffusion for data expressed as a sum of Chebyshev polynomials in 1D
	 * 
	 * This implementation is a simple forward time-difference scheme, 
	 * solving only for the diffusion component of the equations. The class 
	 * solves this implicitly with two matrix equations, one for the even 
	 * and one for the odd Chebyshev polynomials. It uses the BLAS library 
	 * to do so.
	 *********************************************************************/
	class cheb_1D : public plan
	{
	public:
		/*!*******************************************************************
		 * \param i_coeff A double containing the coefficient in front of the diffusion term in the differential equation
		 * \param i_alpha A double that determines the degree of implicit calculation (0.0 = explicit, 1.0 = implicit, 0.5 recommended)
		 * \param i_n An integer number of data elements (grid points) that cheb_1D will be built to tackle
		 * \param i_data_in A double pointer pointing to the input data
		 * \param i_data_out A double pointer pointing to the output data
		 * \param i_flags An integer containing the binary boundary and execution flags
		 *********************************************************************/
		cheb_1D (double i_coeff, double i_alpha, int i_n, double *i_data_in, double *i_data_out = NULL, int flags = 0x00);
		
		virtual ~cheb_1D () {}
		
		/*!*******************************************************************
		 * \brief Execute the operation on the data for a given timestep duration
		 * 
		 * \param timestep a double duration over which the diffusion step will happen
		 *********************************************************************/
		void execute (double timestep);
		
		/*!*******************************************************************
		 * \brief Build either the implicit or explicit matrix for the equation
		 * 
		 * This method takes one of the matrices to be used in calculation and
		 * computes the correct matrix. Since the matrices need only be updated
		 * when the timestep changes, there is potential for a speed increase 
		 * here.
		 * 
		 * \param alpha_scalar The quantity that precedes the diffusion term
		 * \param matrix A double pointer to the matrix to update
		 *********************************************************************/
		void matrix (double alpha_scalar, double *matrix);
		
	private:
		double coeff; //!< A double that represents the coefficient in front of the diffusion term in the differential equation
		double alpha; //!< A double that determines the degree of implicit calculation (0.0 = explicit, 1.0 = implicit, 0.5 recommended)
		int n; //!< An integer number of data elements (grid points) that cheb_1D will be built to handle
		double *data_in; //!< A double pointer to the input data
		double *data_out; //!< A double pointer to the output data; if data_in == data_out, the operation is done in place
		int flags; //!< An integer containing the binary boundary and execution flags
		std::vector<double> diffusion_matrix; //!< A 1D vector to be filled with the implicit matrix equation for the Chebyshev polynomials
		std::vector<double> pre_matrix; //!< A 1D vector to be filled with the explicit matrix equation for the Chebyshev polynomials
		std::vector<double> temp; //!< A temporary 1D vector to store the intermediate step
		std::vector<int> ipiv; //!< An integer vector that contains the reordering for use in the LAPACK routine
		std::unique_ptr<collocation::cheb_grid> cheb; //!< A pointer to a collocation grid that contains the the Chebyshev values
	};
} /* diffusion */

#endif /* end of include guard: diffusion_H_1AIIK1RA */
