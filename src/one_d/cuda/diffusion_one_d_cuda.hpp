/*!**********************************************************************
 * \file diffusion_one_d_cuda.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-26.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DIFFUSION_ONE_D_CUDA_HPP_ORJXCKL7
#define DIFFUSION_ONE_D_CUDA_HPP_ORJXCKL7

#include "../../bases/plan_one_d.hpp"
#include "../../utils/cuda/utils_cublas.hpp"
#include "../../bases/grid.hpp"

namespace cuda
{
	namespace one_d
	{
		/*!*******************************************************************
		* \brief Implementation of explicit diffusion for data expressed as a sum of orthogonal functions in 1D
		* 
		* This implementation calculates the derivatives at the current time
		*********************************************************************/
		template <class datatype>
		class explicit_diffusion : public bases::explicit_plan <datatype>
		{
		public:
			/*!*******************************************************************
			* \param i_coeff A datatype containing the coefficient in front of the diffusion term in the differential equation
			* \param i_grid a shared pointer to the grid, which must be defined for the second derivative
			* \copydoc bases::explicit_plan <datatype>::explicit_plan ()
			*********************************************************************/
			explicit_diffusion (datatype i_coeff, int i_n, bases::grid <datatype>* i_grid, datatype* i_data_in, datatype* i_data_out = NULL);

			virtual ~explicit_diffusion () {}

			/*!*******************************************************************
			* \copydoc bases::explicit_plan <datatype>::execute ()
			*********************************************************************/
			void execute ();

		private:
			using bases::explicit_plan <datatype>::n;
			using bases::explicit_plan <datatype>::data_in;
			using bases::explicit_plan <datatype>::data_out;
		
			datatype coeff; //!< A datatype that represents the coefficient in front of the diffusion term in the differential equation

			utils::vector <datatype> deriv_matrix; //!< A pointer to a collocation grid that contains the the Chebyshev values
		};
	} /* oned */
} /* cuda */

#endif /* end of include guard: DIFFUSION_ONE_D_CUDA_HPP_ORJXCKL7 */
