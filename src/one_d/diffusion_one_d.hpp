/*!***********************************************************************
 * \file diffusion_one_d.hpp
 * Spectral Element
 * 
 * This file provides several implementations of implicit and explicit 
 * methods for solving diffusion.
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DIFFUSION_ONE_D_H_1AIIK1RA
#define DIFFUSION_ONE_D_H_1AIIK1RA

#include <memory>
#include <vector>
#include "../bases/plan.hpp"
#include "../bases/collocation.hpp"
#include "../utils/utils.hpp"

namespace one_d
{
	namespace chebyshev
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
			* \param i_grid a shared pointer to the collocation_grid, which must be defined for the second derivative
			* \copydoc bases::explicit_plan::explicit_plan ()
			*********************************************************************/
			explicit_diffusion (double i_coeff, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, int i_name_in, int i_position, int i_name_out = null);

			virtual ~explicit_diffusion () {}

			/*!*******************************************************************
			* \copydoc bases::explicit_plan::execute ()
			*********************************************************************/
			void execute ();

		private:
			double coeff; //!< A double that represents the coefficient in front of the diffusion term in the differential equation
			int position;
			std::shared_ptr<bases::collocation_grid> grid; //!< A pointer to a collocation grid that contains the the Chebyshev values
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
			 * \param i_timestep A pointer to the double current timestep duration
			 * \param i_grid a shared pointer to the collocation_grid, which must be defined for the second derivative
			 * 
			 * \copydetails bases::implicit_plan::implicit_plan ()
			 *********************************************************************/
			implicit_diffusion (double i_coeff, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_matrix);

			virtual ~implicit_diffusion () {}

			/*!*******************************************************************
			 * \copydoc bases::implicit_plan::execute ()
			 *********************************************************************/
			void execute ();

		private:
			double coeff; //!< A double that represents the coefficient in front of the diffusion term in the differential equation
			std::vector<double> temp;
		};
	} /* chebyshev */
} /* oned */

#endif /* end of include guard: DIFFUSION_ONE_D_H_1AIIK1RA */
