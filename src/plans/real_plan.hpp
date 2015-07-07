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
#include "plan.hpp"

namespace plans
{
	/*!*******************************************************************
	 * \brief A subclass of plan, specific to explicit methods
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
		int dims;
		grids::grid <datatype> &grid_n; //!< A reference to the horizontal grid object
		grids::grid <datatype> &grid_m; //!< A reference to the vertical grid object
		
	public:
		/*!**********************************************************************
		 * \param i_grid_n The grid object associated with the horizontal direction
		 * \param i_grid_m The grid object associated with the vertical direction
		 * \param i_data_in The input data for the plan
		 * \param i_data_out The output data for the plan
		 * \param i_element_flags A pointer to integer flags associated with the element on the whole
		 * \param i_component_flags A pointer to the integer flags associated with the variable associated with the plan
		 ************************************************************************/
		real_plan (grids::variable <datatype> &i_data_in, datatype *i_data_out = NULL, datatype i_coeff = 1.0, int *i_element_flags = NULL, int *i_component_flags = NULL) : 
		plans::plan <datatype> (i_data_in.ptr (), i_data_out, i_element_flags, i_component_flags, i_coeff), 
		n (i_data_in.get_grid (0).get_n ()), 
		ldn (i_data_in.get_grid (0).get_ld ()), 
		m (i_data_in.get_grid (1).get_n ()), 
		dims (i_data_in.dims ()), 
		grid_n (i_data_in.get_grid (0)), 
		grid_m (i_data_in.get_grid (1)) {}

		virtual ~real_plan () {}
		
		virtual void setup () {}
		
		/*!*******************************************************************
		 * \copydoc plans::plan::execute ()
		 *********************************************************************/
		virtual void execute () = 0;
		
		/*!**********************************************************************
		 * \brief An abstract factory class designed to produce a real_plan instance
		 * 
		 * This factory class is implemented to make the creation of plan objects more intuitive. A factory can be given to an equation object, and the equation will know about the grids, data, and flags, so these need not be specified by the user.
		 ************************************************************************/
		class factory : public plan <datatype>::factory
		{
		public:
			factory (datatype i_coeff = 1.0) : plan <datatype>::factory (i_coeff) {}

			virtual ~factory () {}

			virtual const int type () const {
				return plan <datatype>::factory::real;
			}
		};
	};
} /* plans */

#endif /* end of include guard: REAL_PLAN_HPP_CBB7844E */
