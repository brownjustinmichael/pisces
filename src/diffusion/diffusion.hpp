/*!***********************************************************************
 * \file diffusion.hpp
 * Spectral Element
 * 
 * This file provides several implementations of implicit and explicit 
 * methods for solving diffusion.
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
 * \brief Function from BLAS for vector-vector addition (dy = da * dx + dy)
 * 
 * \param n A pointer to the integer number of values in dy
 * \param da A pointer to the double da
 * \param dx The double vector dx
 * \param incx A pointer to the integer spacing of elements in dx
 * \param dy The double vector dy
 * \param incy A pointer to the integer spacing of elements in dy
 *********************************************************************/
extern "C" void daxpy_ (int *n, double *da, double *dx, int *incx, double *dy, int *incy);

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

namespace one_d
{
	namespace explicit_plans
	{
		/*!*******************************************************************
		* \brief Implementation of explicit diffusion for data expressed as a sum of orthogonal functions in 1D
		* 
		* This implementation is a simple forward time-difference scheme, 
		* solving only for the diffusion component of the equations. The class 
		* 
		*********************************************************************/
		class diffusion : public explicit_plan
		{
		public:
			/*!*******************************************************************
			* \param i_coeff A double containing the coefficient in front of the diffusion term in the differential equation
			* \param i_alpha A double that determines the degree of implicit calculation (0.0 = explicit, 1.0 = implicit, 0.5 recommended)
			* \param i_n An integer number of data elements (grid points) that collocation will be built to tackle
			* \param i_data_in A double pointer to the input data
			* \param i_rhs The double array of the right-hand-side, overwritten each timestep with the full right-hand-side (can equal i_data_out)
			* \param i_data_out A double pointer to the output data (if NULL or the same as i_data_in, the operation occurs in place but uses an additional call of dcopy_)
			* \param i_flags_ptr A pointer to an integer containing the binary boundary and execution flags
			*********************************************************************/
			diffusion (double i_coeff, double *i_timestep_ptr, int i_n, std::shared_ptr<collocation::collocation_grid> i_grid, double *i_data_in, double *i_data_out, int *flags_ptr = NULL);
	
			virtual ~diffusion () {}
	
			/*!*******************************************************************
			* @copydoc plan::execute ()
			*********************************************************************/
			void execute ();
		
			/*!*******************************************************************
			* \brief Make a unique pointer to a new collocation object
			* 
			* @copydetails collocation ()
			*********************************************************************/
			inline static std::unique_ptr<plan> make_unique (double i_coeff, double *i_timestep_ptr, int i_n, std::shared_ptr<collocation::collocation_grid> i_grid, double *i_data_in, double *i_data_out, int *flags_ptr = NULL) {
			return std::unique_ptr<plan> (new diffusion (i_coeff, i_timestep_ptr, i_n, i_grid, i_data_in, i_data_out, flags_ptr));
			}
	
		private:
			double coeff; //!< A double that represents the coefficient in front of the diffusion term in the differential equation
			double *timestep_ptr;
			int *flags_ptr; //!< A pointer to an integer containing the binary boundary and execution flags
			std::shared_ptr<collocation::collocation_grid> grid; //!< A pointer to a collocation grid that contains the the Chebyshev values
		};
	} /* explicit_plans */
	
	namespace implicit_plans
	{
		class diffusion : public implicit_plan
		{
		public:
			/*!*******************************************************************
			 * \param i_coeff A double containing the coefficient in front of the diffusion term in the differential equation
			 * \param i_alpha A double that determines the degree of implicit calculation (0.0 = explicit, 1.0 = implicit, 0.5 recommended)
			 * \param i_n An integer number of data elements (grid points) that collocation will be built to tackle
			 * \param i_data_in A double pointer to the input data
			 * \param i_rhs The double array of the right-hand-side, overwritten each timestep with the full right-hand-side (can equal i_data_out)
			 * \param i_data_out A double pointer to the output data (if NULL or the same as i_data_in, the operation occurs in place but uses an additional call of dcopy_)
			 * \param i_flags_ptr A pointer to an integer containing the binary boundary and execution flags
			 *********************************************************************/
			diffusion (double i_coeff, double i_alpha_0, double i_alpha_n, double *i_timestep_ptr, int i_n, std::shared_ptr<collocation::collocation_grid> i_grid, double *i_matrix, int *flags_ptr = NULL);
	
			virtual ~diffusion () {}
	
			/*!*******************************************************************
			 * @copydoc plan::execute ()
			 *********************************************************************/
			void execute ();
		
			/*!*******************************************************************
			 * \brief Make a unique pointer pointing to a new instance of collocation
			 * 
			 * @copydetails collocation ()
			 *********************************************************************/
			inline static std::unique_ptr<plan> make_unique (double i_coeff, double i_alpha_0, double i_alpha_n, double *i_timestep_ptr, int i_n, std::shared_ptr<collocation::collocation_grid> i_grid, double *i_matrix, int *flags_ptr = NULL) {
				return std::unique_ptr<plan> (new diffusion (i_coeff, i_alpha_0, i_alpha_n, i_timestep_ptr, i_n, i_grid, i_matrix, flags_ptr));
			}
	
		private:
			double coeff; //!< A double that represents the coefficient in front of the diffusion term in the differential equation
			double alpha_0;
			double alpha_n;
			double *timestep_ptr;
			int *flags_ptr; //!< A pointer to an integer containing the binary boundary and execution flags
			std::vector<double> temp;
			std::shared_ptr<collocation::collocation_grid> grid; //!< A pointer to a collocation grid that contains the the Chebyshev values
		};
	} /* implicit_plans */
} /* oned */

#endif /* end of include guard: diffusion_H_1AIIK1RA */
