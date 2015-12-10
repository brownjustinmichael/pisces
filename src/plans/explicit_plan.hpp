/*!**********************************************************************
 * \file explicit_plan.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef EXPLICIT_PLAN_HPP_5171CAFD
#define EXPLICIT_PLAN_HPP_5171CAFD

#include "grids/grid.hpp"
#include "plan.hpp"

namespace plans
{
	/*!*******************************************************************
	 * \brief A subclass of plan, specific to explicit methods
	 * 
	 * These plans take input and produce output. Unlike implicit plans, they do not alter the matrix that is solved. Their output also does not need to be retransformed before the solve happens, unlike real plans.
	 *********************************************************************/
	class explicit_plan : public plan
	{
	protected:
		using plan::coeff;
		using plan::data_in;
		using plan::data_out;
		int n; //!< An integer number of data elements (grid points) in the horizontal
		int ldn; //!< The integer max dimension in the horizontal direction
		int m; //!< The integer number of data elements in the vertical
		int dims; //!< The number of dimensions the plan operates on
		grids::grid &grid_n; //!< A reference to the horizontal grid object
		grids::grid &grid_m; //!< A reference to the vertical grid object
		
	public:
		/*!**********************************************************************
		 * \param i_data_in A reference to the input data variable for the plan
		 * \param i_data_out A reference to the output data variable for the plan
		 * @param i_coeff A coefficient by which the plan output should be multiplied
		 ************************************************************************/
		explicit_plan (grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0) : 
		plans::plan (i_data_in, i_data_out, real_spectral, real_spectral, i_coeff), 
		n (i_data_in.get_grid (0).get_n ()), 
		ldn (i_data_in.get_grid (0).get_ld ()), 
		m (i_data_in.get_grid (1).get_n ()), 
		dims (i_data_in.dims ()), 
		grid_n (i_data_in.get_grid (0)), 
		grid_m (i_data_in.get_grid (1)) {}

		virtual ~explicit_plan () {}
				
		/*!*******************************************************************
		 * \copydoc plans::plan::execute ()
		 *********************************************************************/
		virtual void execute () = 0;
	};
} /* plans */

#endif /* end of include guard: EXPLICIT_PLAN_HPP_5171CAFD */
