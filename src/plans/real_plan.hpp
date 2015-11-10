/*!**********************************************************************
 * \file real_plan.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef REAL_PLAN_HPP_CBB7844E
#define REAL_PLAN_HPP_CBB7844E

#include "grids/grid.hpp"
#include "linalg/utils.hpp"
#include "plan.hpp"

namespace plans
{
	/*!*******************************************************************
	 * \brief A subclass of plan, specific to real methods
	 * 
	 * These plans take input and produce output. Unlike implicit plans, they do not make any changes to the matrices associated with the solve. However, they do require a transform before the solve to put them in the correct state. This transform is not done within the class but rather in equation right before the solve happens (this is to make sure that the transform is done only once).
	 *********************************************************************/
	template <class datatype>
	class real_plan : public plan <datatype>
	{
	protected:
		using plan <datatype>::coeff;
		using plan <datatype>::data_in;
		using plan <datatype>::data_out;
		int n; //!< An integer number of data elements (grid points) in the horizontal
		int ldn; //!< The integer max dimension in the horizontal direction
		int m; //!< The integer number of data elements in the vertical
		int dims; //!< The number of dimensions the plan operates in
		grids::grid &grid_n; //!< A reference to the horizontal grid object
		grids::grid &grid_m; //!< A reference to the vertical grid object
		
	public:
		/*!**********************************************************************
		 * \param i_data_in A reference to the input data variable for the plan
		 * \param i_data_out A reference to the output data variable for the plan
		 * @param i_coeff A coefficient by which the plan output should be multiplied
		 ************************************************************************/
		real_plan (grids::variable &i_data_in, grids::variable &i_data_out, datatype i_coeff = 1.0) : 
		plans::plan <datatype> (i_data_in, i_data_out, real_real, real_real, i_coeff), 
		n (i_data_in.get_grid (0).get_n ()), 
		ldn (i_data_in.get_grid (0).get_ld ()), 
		m (i_data_in.get_grid (1).get_n ()), 
		dims (i_data_in.dims ()), 
		grid_n (i_data_in.get_grid (0)), 
		grid_m (i_data_in.get_grid (1)) {}

		virtual ~real_plan () {}
		
		/*!*******************************************************************
		 * \copydoc plans::plan::execute ()
		 *********************************************************************/
		virtual void execute () = 0;
	};
} /* plans */

#endif /* end of include guard: REAL_PLAN_HPP_CBB7844E */
