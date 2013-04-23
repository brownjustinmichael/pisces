/*!***********************************************************************
 * \file one_d/diffusion.hpp
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
#include "../bases/plan.hpp"
#include "../bases/collocation.hpp"

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
	/*!*******************************************************************
	* \brief Implementation of explicit diffusion for data expressed as a sum of orthogonal functions in 1D
	* 
	* This implementation calculates the derivatives at the current time
	*********************************************************************/
	class explicit_diffusion : public bases::explicit_plan
	{
	public:
		/*!*******************************************************************
		* \param i_coeff A double containing the coefficient in front of the diffusion term in the differential equation
		* \param i_timestep_ptr A pointer to the double current timestep duration
		* \param i_grid a shared pointer to the collocation_grid, which must be defined for the second derivative
		* 
		* \copydetails bases::explicit_plan::explicit_plan ()
		*********************************************************************/
		explicit_diffusion (double i_coeff, double *i_timestep_ptr, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_data_in, double *i_data_out = NULL, int *i_flags_ptr = NULL, int i_logger = -1) : bases::explicit_plan (i_n, i_data_in, i_data_out, i_flags_ptr, i_logger) {
			init (i_coeff, i_timestep_ptr, i_n, i_grid, i_data_in, i_data_out, i_flags_ptr);
		}
		
		explicit_diffusion (double i_coeff, double *i_timestep_ptr, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double &i_data_in, double &i_data_out, int *i_flags_ptr = NULL, int i_logger = -1) : bases::explicit_plan (i_n, &i_data_in, &i_data_out, i_flags_ptr, i_logger) {
			init (i_coeff, i_timestep_ptr, i_n, i_grid, &i_data_in, &i_data_out, i_flags_ptr);
		}
		
		explicit_diffusion (double i_coeff, double *i_timestep_ptr, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double &i_data_in, int *i_flags_ptr = NULL, int i_logger = -1) : bases::explicit_plan (i_n, &i_data_in, NULL, i_flags_ptr, i_logger) {
			init (i_coeff, i_timestep_ptr, i_n, i_grid, &i_data_in);
		}

		virtual ~explicit_diffusion () {}

		/*!*******************************************************************
		* \copydoc bases::explicit_plan::execute ()
		*********************************************************************/
		void execute ();
	
		/*!*******************************************************************
		* \brief Make a unique pointer to a new collocation object
		* 
		* \copydetails explicit_diffusion ()
		*********************************************************************/
		inline static std::unique_ptr<plan> make_unique (double i_coeff, double *i_timestep_ptr, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_data_in, double *i_data_out, int *i_flags_ptr = NULL, int i_logger = -1) {
		return std::unique_ptr<plan> (new explicit_diffusion (i_coeff, i_timestep_ptr, i_n, i_grid, i_data_in, i_data_out, i_flags_ptr, i_logger));
		}
		
		inline static std::unique_ptr<plan> make_unique (double i_coeff, double *i_timestep_ptr, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double &i_data_in, double &i_data_out, int *i_flags_ptr = NULL, int i_logger = -1) {
		return std::unique_ptr<plan> (new explicit_diffusion (i_coeff, i_timestep_ptr, i_n, i_grid, i_data_in, i_data_out, i_flags_ptr, i_logger));
		}
		
		inline static std::unique_ptr<plan> make_unique (double i_coeff, double *i_timestep_ptr, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double &i_data_in, int *i_flags_ptr = NULL, int i_logger = -1) {
		return std::unique_ptr<plan> (new explicit_diffusion (i_coeff, i_timestep_ptr, i_n, i_grid, i_data_in, i_flags_ptr, i_logger));
		}

	private:
		double coeff; //!< A double that represents the coefficient in front of the diffusion term in the differential equation
		double *timestep_ptr;
		int *flags_ptr; //!< A pointer to an integer containing the binary boundary and execution flags
		std::shared_ptr<bases::collocation_grid> grid; //!< A pointer to a collocation grid that contains the the Chebyshev values
		
		void init (double i_coeff, double *i_timestep_ptr, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_data_in, double *i_data_out = NULL, int *i_flags_ptr = NULL);
	};

	/*!*******************************************************************
	 * \brief Implementation of implicit diffusion for data expressed as a sum of orthogonal functions in 1D
	 * 
	 * This implementation adds the diffusion terms to the implicit matrix.
	 *********************************************************************/
	class implicit_diffusion : public bases::implicit_plan
	{
	public:
		/*!*******************************************************************
		 * \param i_coeff A double containing the coefficient in front of the diffusion term in the differential equation
		 * \param i_alpha_0 A double that determines the multiplier on the top boundary
		 * \param i_alpha_n A double that determines the multiplier on the bottom boundary
		 * \param i_timestep_ptr A pointer to the double current timestep duration
		 * \param i_grid a shared pointer to the collocation_grid, which must be defined for the second derivative
		 * 
		 * \copydetails bases::implicit_plan::implicit_plan ()
		 *********************************************************************/
		implicit_diffusion (double i_coeff, double i_alpha_0, double i_alpha_n, double *i_timestep_ptr, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_matrix, int *i_flags_ptr = NULL, int i_logger = -1);

		virtual ~implicit_diffusion () {}

		/*!*******************************************************************
		 * \copydoc bases::implicit_plan::execute ()
		 *********************************************************************/
		void execute ();
	
		/*!*******************************************************************
		 * \brief Make a unique pointer pointing to a new instance of collocation
		 * 
		 * \copydetails implicit_diffusion ()
		 *********************************************************************/
		inline static std::unique_ptr<plan> make_unique (double i_coeff, double i_alpha_0, double i_alpha_n, double *i_timestep_ptr, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_matrix, int *i_flags_ptr = NULL, int i_logger = -1) {
			return std::unique_ptr<plan> (new implicit_diffusion (i_coeff, i_alpha_0, i_alpha_n, i_timestep_ptr, i_n, i_grid, i_matrix, i_flags_ptr, i_logger));
		}

	private:
		double coeff; //!< A double that represents the coefficient in front of the diffusion term in the differential equation
		double alpha_0;
		double alpha_n;
		double *timestep_ptr;
		int *flags_ptr; //!< A pointer to an integer containing the binary boundary and execution flags
		std::vector<double> temp;
		std::shared_ptr<bases::collocation_grid> grid; //!< A pointer to a collocation grid that contains the the Chebyshev values
	};
} /* oned */

#endif /* end of include guard: diffusion_H_1AIIK1RA */
