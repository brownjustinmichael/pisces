/*!**********************************************************************
 * \file plan_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-11.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef IMPLICIT_PLAN_TWO_D_HPP_XQ7AJI7K
#define IMPLICIT_PLAN_TWO_D_HPP_XQ7AJI7K

#include "grid.hpp"
#include "plan.hpp"

namespace plans
{
	namespace two_d
	{
		/*!*******************************************************************
		 * \brief A subclass of plan, specific to implicit methods
		 * 
		 * These plans produce output in a square matrix but take no input.
		 *********************************************************************/
		template <class datatype>
		class implicit_plan : public plans::implicit_plan <datatype>
		{
		public:
			implicit_plan (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, datatype *i_matrix_n, datatype *i_matrix_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) :
			plans::implicit_plan <datatype> (i_element_flags, i_component_flags),  
			matrix_n (i_matrix_n),  
			matrix_m (i_matrix_m),
			n (i_grid_n.get_n ()),
			ldn (i_grid_n.get_ld ()),
			m (i_grid_m.get_n ()),
			grid_n (i_grid_n),
			grid_m (i_grid_m),
			data_in (i_data_in),
			data_out (i_data_out ? i_data_out : i_data_in) {}
		
			virtual ~implicit_plan () {}

			/*!*******************************************************************
			 * \copydoc plan::execute ()
			 *********************************************************************/
			virtual void execute () = 0;
		
			using plans::implicit_plan <datatype>::factory;

			datatype *matrix_n; //!< A datatype pointer to the input data
			datatype *matrix_m; //!< A datatype pointer to the input data
			int n; //!< An integer number of data elements (grid points) that collocation_1D will be built to handle
			int ldn;
			int m;
			plans::grid <datatype> &grid_n, &grid_m;
			datatype* data_in; //!< A datatype pointer to the input data
			datatype* data_out; //!< A datatype pointer to the output data
		};
	} /* two_d */
} /* plans */

#endif /* end of include guard: IMPLICIT_PLAN_TWO_D_HPP_XQ7AJI7K */
