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
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param dx The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 * \param dy The array to which the data are copied
 * \param incy A pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void dcopy_(int *n, double *x, int *incx, double *y, int *incy);

/*!*******************************************************************
 * \brief Function from BLAS that swaps a double array with another in place
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param x The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 * \param y The array to which the data are copied
 * \param incy A pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void dswap_(int *n, double *dx, int *incx, double *dy, int *incy);

extern "C" void daxpy_ (int *n, double *da, double *dx, int *incx, double *y, int *incy);

extern "C" void dscal_ (int *n, double *da, double *dx, int *incx);

/*!*******************************************************************
 * \brief Function from BLAS for matrix-vector multiplication (y = alpha * a * x + beta * y)
 * 
 * \param trans A pointer to transposition character ("N" for not transposed, "T" for transposed)
 * \param m A pointer to the number of rows in a
 * \param n A pointer to the number of columns in a
 * \param alpha A pointer to the double multiplier on a
 * \param a The double matrix a
 * \param lda A pointer to the integer number of leading dimension of a
 * \param x The double vector x
 * \param incx A pointer to an integer spacing of elements in x
 * \param beta A pointer to the double multiplier on y
 * \param y The double vector y, overwritten with the solution
 * \param incy A pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void dgemv_ (char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

/*!*******************************************************************
 * \brief Function from LAPACK that factorizes the matrix a by LU decomposition
 * 
 * \param m A pointer to the number of rows in a
 * \param n A pointer to the number of columns in a
 * \param a A double matrix to be overwritten with its LU decomposition
 * \param ipiv An integer array to contain the pivot indices
 * \param info A pointer to an integer indicating success (0 for successful exit)
 *********************************************************************/
extern "C" void dgetrf_ (int *m, int *n, double *a, int *lda, int *ipiv, int *info);

/*!*******************************************************************
 * \brief Function from LAPACK that solves a factorized matrix equation
 * 
 * \param trans A pointer to transposition character ("N" for not transposed, "T" for transposed)
 * \param n A pointer to the number of columns in a
 * \param nrhs A pointer to the number of right hand sides
 * \param a A double matrix to be overwritten with its LU decomposition
 * \param lda A pointer to the integer number of leading dimension of a
 * \param ipiv An integer array to contain the pivot indices
 * \param b The double right hand side array, overwritten with solution
 * \param ldb A pointer to the integer number of leading dimension of b
 * \param info A pointer to an integer indicating success (0 for successful exit)
 *********************************************************************/
extern "C" void dgetrs_ (char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

extern "C" void dgesv_ (int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

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
	class collocation_chebyshev_explicit_1D : public plan
	{
	public:
		/*!*******************************************************************
		 * \param i_coeff A double containing the coefficient in front of the diffusion term in the differential equation
		 * \param i_alpha A double that determines the degree of implicit calculation (0.0 = explicit, 1.0 = implicit, 0.5 recommended)
		 * \param i_n An integer number of data elements (grid points) that collocation_chebyshev_1D will be built to tackle
		 * \param i_data_in A double pointer to the input data
		 * \param i_rhs The double array of the right-hand-side, overwritten each timestep with the full right-hand-side (can equal i_data_out)
		 * \param i_data_out A double pointer to the output data (if NULL or the same as i_data_in, the operation occurs in place but uses an additional call of dcopy_)
		 * \param i_flags An integer containing the binary boundary and execution flags
		 *********************************************************************/
		collocation_chebyshev_explicit_1D (double i_coeff, int i_n, std::shared_ptr<collocation::chebyshev_grid> i_cheb, double *i_data_in, double *i_data_out, int flags = 0x00);
		
		virtual ~collocation_chebyshev_explicit_1D () {}
		
		/*!*******************************************************************
		 * \brief Execute the operation on the data for a given timestep duration
		 * 
		 * Of particular note, rhs is overwritten with the full right-hand-side 
		 * (including diffusion terms) after execution.
		 * 
		 * \param timestep a double duration over which the diffusion step will happen
		 *********************************************************************/
		void execute (double timestep);
		
	private:
		double coeff; //!< A double that represents the coefficient in front of the diffusion term in the differential equation
		double previous_timestep; //!< A double that records the previous timestep 
		int n; //!< An integer number of data elements (grid points) that collocation_chebyshev_1D will be built to handle
		double *data_in; //!< A double pointer to the input data
		double *data_out; //!< A double pointer to the output data; if data_in == data_out, the operation is done in place (but inefficient)
		int flags; //!< An integer containing the binary boundary and execution flags
		std::shared_ptr<collocation::chebyshev_grid> cheb; //!< A pointer to a collocation grid that contains the the Chebyshev values
	};
	
	class collocation_chebyshev_implicit_1D : public plan
	{
	public:
		/*!*******************************************************************
		 * \param i_coeff A double containing the coefficient in front of the diffusion term in the differential equation
		 * \param i_alpha A double that determines the degree of implicit calculation (0.0 = explicit, 1.0 = implicit, 0.5 recommended)
		 * \param i_n An integer number of data elements (grid points) that collocation_chebyshev_1D will be built to tackle
		 * \param i_data_in A double pointer to the input data
		 * \param i_rhs The double array of the right-hand-side, overwritten each timestep with the full right-hand-side (can equal i_data_out)
		 * \param i_data_out A double pointer to the output data (if NULL or the same as i_data_in, the operation occurs in place but uses an additional call of dcopy_)
		 * \param i_flags An integer containing the binary boundary and execution flags
		 *********************************************************************/
		collocation_chebyshev_implicit_1D (double i_coeff, double i_alpha, int i_n, std::shared_ptr<collocation::chebyshev_grid> i_cheb, double *i_data_in, double *i_rhs, double *i_data_out, int flags = 0x00);
		
		virtual ~collocation_chebyshev_implicit_1D () {}
		
		/*!*******************************************************************
		 * \brief Execute the operation on the data for a given timestep duration
		 * 
		 * Of particular note, rhs is overwritten with the full right-hand-side 
		 * (including diffusion terms) after execution.
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
		double previous_timestep; //!< A double that records the previous timestep 
		int n; //!< An integer number of data elements (grid points) that collocation_chebyshev_1D will be built to handle
		double *data_in; //!< A double pointer to the input data
		double *data_out; //!< A double pointer to the output data; if data_in == data_out, the operation is done in place (but inefficient)
		double *rhs; //!< A double pointer to the right-hand-side; if rhs == data_out, the operation is most efficient
		int flags; //!< An integer containing the binary boundary and execution flags
		std::vector<double> temp;
		std::vector<double> diffusion_matrix; //!< A 1D vector to be filled with the implicit matrix equation for the Chebyshev polynomials
		std::vector<double> pre_matrix; //!< A 1D vector to be filled with the explicit matrix equation for the Chebyshev polynomials
		std::vector<int> ipiv; //!< An integer vector that contains the reordering for use in the LAPACK routine
		std::shared_ptr<collocation::chebyshev_grid> cheb; //!< A pointer to a collocation grid that contains the the Chebyshev values
	};
} /* diffusion */

#endif /* end of include guard: diffusion_H_1AIIK1RA */
