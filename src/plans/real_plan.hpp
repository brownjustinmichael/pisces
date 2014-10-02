/*!**********************************************************************
 * \file real_plan.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef REAL_PLAN_HPP_CBB7844E
#define REAL_PLAN_HPP_CBB7844E

#include "grid.hpp"
#include "plan.hpp"

namespace plans
{
	namespace two_d
	{
		/*!*******************************************************************
		 * \brief A subclass of plan, specific to explicit methods
		 * 
		 * These plans take input and produce output.
		 *********************************************************************/
		template <class datatype>
		class real_plan : public plans::real_plan <datatype>
		{
		public:
			real_plan (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) :
			plans::real_plan <datatype> (i_element_flags, i_component_flags),
			n (i_grid_n.get_n ()),
			ldn (i_grid_n.get_ld ()),
			m (i_grid_m.get_n ()),
			grid_n (i_grid_n),
			grid_m (i_grid_m),
			data_in (i_data_in),
			data_out (i_data_out ? i_data_out : i_data_in) {}

			virtual ~real_plan () {}

			/*!*******************************************************************
			 * \copydoc plans::plan::execute ()
			 *********************************************************************/
			virtual void execute () = 0;
	
			int n; //!< An integer number of data elements (grid points) that collocation_1D will be built to handle
			int ldn;
			int m;
			plans::grid <datatype> &grid_n, &grid_m;
			datatype* data_in; //!< A datatype pointer to the input data
			datatype* data_out; //!< A datatype pointer to the output data
		};
	} /* two_d */
} /* plans */

#endif /* end of include guard: REAL_PLAN_HPP_CBB7844E */
