/*!**********************************************************************
 * \file implicit_plan.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-11.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef IMPLICIT_PLAN_TWO_D_HPP_XQ7AJI7K
#define IMPLICIT_PLAN_TWO_D_HPP_XQ7AJI7K

#include "grids/grid.hpp"
#include "plan.hpp"

namespace plans
{
	/*!*******************************************************************
	 * \brief A subclass of plan, specific to implicit methods
	 * 
	 * These plans take input and give output, just like real and explicit plans; however, they also have access to the matrices associated with the solves. This lets any subclasses add values to those matrices. Any matrix editing should happen in the constructor. The execute method here should be reserved for plans that need to have some component evaluated explicitly, such as in the case of alternating direction implicit solves, where every half timestep, the matrix is evaluated explicitly.
	 *********************************************************************/
	template <class datatype>
	class implicit_plan : public plan <datatype>
	{
	protected:
		datatype *matrix_n; //!< A datatype pointer to the input data
		datatype *matrix_m; //!< A datatype pointer to the input data
		int n; //!< An integer number of data elements (grid points) that collocation_1D will be built to handle
		int ldn; //!< The integer max dimension in the horizontal direction
		int m; //!< The integer number of data elements in the vertical
		grids::grid <datatype> &grid_n; //!< A reference to the horizontal grid object
		grids::grid <datatype> &grid_m; //!< A reference to the vertical grid object
		datatype* data_in; //!< A datatype pointer to the input data
		datatype* data_out; //!< A datatype pointer to the output data
		
	public:
		/*!**********************************************************************
		 * \param i_grid_n The grid object associated with the horizontal direction
		 * \param i_grid_m The grid object associated with the vertical direction
		 * \param i_matrix_n The datatype matrix associated with the horizontal solve
		 * \param i_matrix_m The datatype matrix associated with the vertical solve
		 * \param i_data_in The input data for the plan
		 * \param i_data_out The output data for the plan
		 * \param i_element_flags A pointer to integer flags associated with the element on the whole
		 * \param i_component_flags A pointer to the integer flags associated with the variable associated with the plan
		 ************************************************************************/
		implicit_plan (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype *i_matrix_n, datatype *i_matrix_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) :
		plans::plan <datatype> (i_element_flags, i_component_flags), matrix_n (i_matrix_n), matrix_m (i_matrix_m), n (i_grid_n.get_n ()), ldn (i_grid_n.get_ld ()), m (i_grid_m.get_n ()), grid_n (i_grid_n), grid_m (i_grid_m), data_in (i_data_in), data_out (i_data_out ? i_data_out : i_data_in) {}
	
		virtual ~implicit_plan () {}

		/*!*******************************************************************
		 * \copydoc plan::execute ()
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
			 * \param matrices An array of the solver datatype matrices
			 * \param i_data_in The datatype array of the input data
			 * \param i_data_out The datatype array of the output data
			 * \param i_element_flags A pointer to the integer flags associated with the element on the whole
			 * \param i_component_flags A pointer to the integer flags associated with the variable associated with the plan
			 * 
			 * This method creates a shared_ptr to an implicit plan instance. The benefit to this inclusion is that the instance method can be called in a uniform way and hide communication of grid and matrix information from the user. If a plan would be created that would not do anything (e.g. something with a coefficient of 0.0), this will return a NULL shared pointer.
			 ************************************************************************/
			virtual std::shared_ptr <plans::plan <datatype>> instance (grids::grid <datatype> **grids, datatype **matrices, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const = 0;
		};
	};
} /* plans */

#endif /* end of include guard: IMPLICIT_PLAN_TWO_D_HPP_XQ7AJI7K */
