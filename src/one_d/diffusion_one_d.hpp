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
			* \param i_timestep A pointer to the double current timestep duration
			* \param i_grid a shared pointer to the collocation_grid, which must be defined for the second derivative
			* 
			* \copydetails bases::explicit_plan::explicit_plan ()
			*********************************************************************/
			explicit_diffusion (double i_coeff, double& i_timestep, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_data_in, double *i_data_out = NULL, int *i_flags_ptr = NULL, int i_logger = -1) : bases::explicit_plan (i_n, i_data_in, i_data_out, i_flags_ptr, i_logger), timestep (i_timestep) {
				init (i_coeff, i_timestep, i_n, i_grid, i_data_in, i_data_out, i_flags_ptr);
			}
		
			explicit_diffusion (double i_coeff, double& i_timestep, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double& i_data_in, double& i_data_out, int *i_flags_ptr = NULL, int i_logger = -1) : bases::explicit_plan (i_n, &i_data_in, &i_data_out, i_flags_ptr, i_logger), timestep (i_timestep) {
				init (i_coeff, i_timestep, i_n, i_grid, &i_data_in, &i_data_out, i_flags_ptr);
			}
		
			explicit_diffusion (double i_coeff, double& i_timestep, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double& i_data_in, int *i_flags_ptr = NULL, int i_logger = -1) : bases::explicit_plan (i_n, &i_data_in, NULL, i_flags_ptr, i_logger), timestep (i_timestep) {
				init (i_coeff, i_timestep, i_n, i_grid, &i_data_in);
			}

			virtual ~explicit_diffusion () {}

			/*!*******************************************************************
			* \copydoc bases::explicit_plan::execute ()
			*********************************************************************/
			void execute ();

		private:
			double coeff; //!< A double that represents the coefficient in front of the diffusion term in the differential equation
			double& timestep;
			int *flags_ptr; //!< A pointer to an integer containing the binary boundary and execution flags
			std::shared_ptr<bases::collocation_grid> grid; //!< A pointer to a collocation grid that contains the the Chebyshev values
		
			void init (double i_coeff, double& i_timestep, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_data_in, double *i_data_out = NULL, int *i_flags_ptr = NULL);
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
			implicit_diffusion (double i_coeff, double& i_timestep, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_matrix, int *i_flags_ptr = NULL, int i_logger = -1);

			virtual ~implicit_diffusion () {}

			/*!*******************************************************************
			 * \copydoc bases::implicit_plan::execute ()
			 *********************************************************************/
			void execute ();

		private:
			double coeff; //!< A double that represents the coefficient in front of the diffusion term in the differential equation
			double& timestep;
			int *flags_ptr; //!< A pointer to an integer containing the binary boundary and execution flags
			std::vector<double> temp;
		};
	} /* chebyshev */
} /* oned */

#endif /* end of include guard: DIFFUSION_ONE_D_H_1AIIK1RA */
