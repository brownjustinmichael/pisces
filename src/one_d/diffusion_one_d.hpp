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
		* \param i_grid a shared pointer to the collocation_grid, which must be defined for the second derivative
		* \copydoc bases::explicit_plan <datatype>::explicit_plan ()
		*********************************************************************/
		explicit_diffusion (bases::element <datatype>* i_element_ptr, datatype i_coeff, int i_n, bases::collocation_grid <datatype>* i_grid, datatype* i_data_in, datatype* i_data_out = null, int i_flags = 0x00);

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
		bases::collocation_grid <datatype>* grid; //!< A pointer to a collocation grid that contains the the Chebyshev values
	};

	/*!*******************************************************************
	 * \brief Implementation of implicit diffusion for data expressed as a sum of orthogonal functions in 1D
	 * 
	 * This implementation adds the diffusion terms to the implicit matrix.
	 *********************************************************************/
	template <class datatype>
	class implicit_diffusion : public bases::implicit_plan <datatype>
	{
	public:
		/*!*******************************************************************
		 * \param i_coeff A datatype containing the coefficient in front of the diffusion term in the differential equation
		 * \param i_grid a shared pointer to the collocation_grid, which must be defined for the second derivative
		 * 
		 * \copydetails bases::implicit_plan <datatype>::implicit_plan ()
		 *********************************************************************/
		implicit_diffusion (bases::element <datatype>* i_element_ptr, datatype i_coeff, int i_n, bases::collocation_grid <datatype>* i_grid, datatype *i_matrix, int i_flags = 0x00);

		virtual ~implicit_diffusion () {}

		/*!*******************************************************************
		 * \copydoc bases::implicit_plan <datatype>::execute ()
		 *********************************************************************/
		void execute ();

	private:
		using bases::implicit_plan <datatype>::n;
		using bases::implicit_plan <datatype>::grid;
		using bases::implicit_plan <datatype>::matrix;
		
		datatype coeff; //!< A datatype that represents the coefficient in front of the diffusion term in the differential equation
	};
} /* oned */

#endif /* end of include guard: DIFFUSION_ONE_D_H_1AIIK1RA */
