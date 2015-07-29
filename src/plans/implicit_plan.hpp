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
#include "io/parameters.hpp"

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
		using plan <datatype>::coeff;
		using plan <datatype>::data_in;
		using plan <datatype>::data_out;
		datatype *matrix_n; //!< A datatype pointer to the input data
		datatype *matrix_m; //!< A datatype pointer to the input data
		int n; //!< An integer number of data elements (grid points) that collocation_1D will be built to handle
		int ldn; //!< The integer max dimension in the horizontal direction
		int m; //!< The integer number of data elements in the vertical
		int dims;
		grids::grid <datatype> &grid_n; //!< A reference to the horizontal grid object
		grids::grid <datatype> &grid_m; //!< A reference to the vertical grid object
		
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
		implicit_plan (datatype *i_matrix_n, datatype *i_matrix_m, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0, int state_in = spectral_spectral, int state_out = spectral_spectral) :
		plans::plan <datatype> (i_data_in, i_data_out, state_in, state_out, i_coeff), 
		matrix_n (i_matrix_n), 
		matrix_m (i_matrix_m), 
		n (i_data_in.get_grid (0).get_n ()), 
		ldn (i_data_in.get_grid (0).get_ld ()), 
		m (i_data_in.get_grid (1).get_n ()), 
		dims (i_data_in.dims ()), 
		grid_n (i_data_in.get_grid (0)), 
		grid_m (i_data_in.get_grid (1)) {}
	
		virtual ~implicit_plan () {}

		virtual bool implicit () {
			return true;
		}
		
		virtual void setup () = 0;
		
		/*!*******************************************************************
		 * \copydoc plan::execute ()
		 *********************************************************************/
		virtual void execute () = 0;

		/*!**********************************************************************
		 * \brief An abstract factory class designed to produce an implicit_plan instance
		 * 
		 * This factory class is implemented to make the creation of plan objects more intuitive. A factory can be given to an equation object, and the equation will know about the grids, data, and flags, so these need not be specified by the user.
		 ************************************************************************/
		class factory : public plan <datatype>::factory
		{
		public:
			factory (datatype i_coeff = 1.0) : plan <datatype>::factory (i_coeff) {}

			virtual ~factory () {}
		};
	};
} /* plans */

#endif /* end of include guard: IMPLICIT_PLAN_TWO_D_HPP_XQ7AJI7K */
