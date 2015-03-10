/*!**********************************************************************
 * \file explicit_plan.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef EXPLICIT_PLAN_HPP_5171CAFD
#define EXPLICIT_PLAN_HPP_5171CAFD

#include "grid.hpp"
#include "plan.hpp"

namespace plans
{
	/*!*******************************************************************
	 * \brief A subclass of plan, specific to explicit methods
	 * 
	 * These plans take input and produce output. Unlike implicit plans, they do not alter the matrix that is solved. Their output also does not need to be retransformed before the solve happens, unlike real plans.
	 *********************************************************************/
	template <class datatype>
	class explicit_plan : public plan <datatype>
	{
	protected:
		int n; //!< An integer number of data elements (grid points) in the horizontal
		int ldn; //!< The integer max dimension in the horizontal direction
		int m; //!< The integer number of data elements in the vertical
		plans::grid <datatype> &grid_n; //!< A reference to the horizontal grid object
		plans::grid <datatype> &grid_m; //!< A reference to the vertical grid object
		datatype* data_in; //!< A datatype pointer to the input data
		datatype* data_out; //!< A datatype pointer to the output data
		
	public:
		/*!**********************************************************************
		 * \param i_grid_n The grid object associated with the horizontal direction
		 * \param i_grid_m The grid object associated with the vertical direction
		 * \param i_data_in The input data for the plan
		 * \param i_data_out The output data for the plan
		 * \param i_element_flags A pointer to integer flags associated with the element on the whole
		 * \param i_component_flags A pointer to the integer flags associated with the variable associated with the plan
		 ************************************************************************/
		explicit_plan (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) : plans::plan <datatype> (i_element_flags, i_component_flags), n (i_grid_n.get_n ()), ldn (i_grid_n.get_ld ()), m (i_grid_m.get_n ()), grid_n (i_grid_n), grid_m (i_grid_m), data_in (i_data_in), data_out (i_data_out ? i_data_out : i_data_in) {}

		virtual ~explicit_plan () {}

		/*!*******************************************************************
		 * \copydoc plans::plan::execute ()
		 *********************************************************************/
		virtual void execute () = 0;
	
		/*!**********************************************************************
		 * \brief An abstract factory class designed to produce an explicit_plan instance
		 * 
		 * This factory class is implemented to make the creation of plan objects more intuitive. A factory can be given to an equation object, and the equation will know about the grids, data, and flags, so these need not be specified by the user.
		 ************************************************************************/
		class factory
		{
		public:
			virtual ~factory () {}
			
			/*!**********************************************************************
			 * \brief The abstract instance creating class
			 * 
			 * \param grids An array of grid objects that define the data grid
			 * \param i_data_in The datatype array of the input data
			 * \param i_data_out The datatype array of the output data
			 * \param i_element_flags A pointer to the integer flags associated with the element on the whole
			 * \param i_component_flags A pointer to the integer flags associated with the variable associated with the plan
			 * 
			 * This method creates a shared_ptr to an explicit plan instance. The benefit to this inclusion is that the instance method can be called in a uniform way and hide communication of grid and matrix information from the user. If a plan would be created that would not do anything (e.g. something with a coefficient of 0.0), this will return a NULL shared pointer.
			 ************************************************************************/
			virtual std::shared_ptr <plans::plan <datatype>> instance (plans::grid <datatype> **grids, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const = 0;
		};
	};
} /* plans */

#endif /* end of include guard: EXPLICIT_PLAN_HPP_5171CAFD */
