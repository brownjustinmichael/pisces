/*!**********************************************************************
 * \file plan_one_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-11.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef PLAN_HPP_8B98W8AH
#define PLAN_HPP_8B98W8AH

#include "grid.hpp"
#include "plan.hpp"
#include "../bases/solver.hpp"

namespace one_d
{
	namespace chebyshev
	{
		enum mode {
			mode_flag = 0x01
		};
	} /* chebyshev */
	
	namespace cosine
	{
		enum mode {
			mode_flag = 0x02
		};
	} /* chebyshev */
	
	/*!*******************************************************************
	 * \brief A subclass of plan, specific to explicit methods
	 * 
	 * These plans take input and produce output.
	 *********************************************************************/
	template <class datatype>
	class explicit_plan : public bases::plan <datatype>
	{
	public:
		explicit_plan (bases::grid <datatype> &i_grid, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) : 
		bases::plan <datatype> (i_element_flags, i_component_flags),
		n (i_grid.get_n ()),
		ld (i_grid.get_ld ()),
		grid (i_grid),
		data_in (i_data_in),
		data_out (i_data_out ? i_data_out : i_data_in) {}
		
		/*!*******************************************************************
		 * \param i_n The integer number of elements in the data
		 * \param i_data_in The integer scalar index of the input
		 * \param i_data_out The integer scalar index of the output
		 * \copydoc plan::plan ()
		 *********************************************************************/
		explicit_plan (bases::master_solver <datatype> &i_solver) : 
		bases::plan <datatype> (i_solver.element_flags, i_solver.component_flags),
		n (i_solver.grid_ptr ()->get_n ()),
		ld (i_solver.grid_ptr ()->get_ld ()),
		grid (*(i_solver.grid_ptr ())),
		data_in (i_solver.data_ptr ()),
		data_out (i_solver.rhs_ptr (explicit_rhs)) {}

		virtual ~explicit_plan () {
			// printf ("Destroying one_d explicit plan\n");
		}

		/*!*******************************************************************
		 * \copydoc bases::plan::execute ()
		 *********************************************************************/
		virtual void execute () = 0;

	protected:

		int n, ld; //!< An integer number of data elements (grid points) that collocation_1D will be built to handle
		bases::grid <datatype> &grid;
		datatype* data_in; //!< A datatype pointer to the input data
		datatype* data_out; //!< A datatype pointer to the output data
	};

	/*!*******************************************************************
	 * \brief A subclass of plan, specific to implicit methods
	 * 
	 * These plans produce output in a square matrix but take no input.
	 *********************************************************************/
	template <class datatype>
	class implicit_plan : public explicit_plan <datatype>
	{
	public:
		implicit_plan (bases::grid <datatype> &i_grid, datatype *i_matrix, datatype *i_data_in, datatype *i_data_out = NULL, int *element_flags = NULL, int *component_flags = NULL) :
		explicit_plan <datatype> (i_grid, i_data_in, i_data_out),  
		matrix (i_matrix) {}
		
		/*!*******************************************************************
		 * \param i_n The integer number of elements in a row of the square i_matrix
		 * \param i_grid A shared pointer to the collocation grid object
		 * \param i_matrix The datatype matrix to be updated
		 * \copydoc plan::plan ()
		 *********************************************************************/
		implicit_plan (bases::master_solver <datatype> &i_solver) :
		explicit_plan <datatype> (i_solver),  
		matrix (i_solver.matrix_ptr ()) {
			data_out = i_solver.rhs_ptr (0);
		}

		virtual ~implicit_plan () {
			// printf ("Destroying one_d implicit plan\n");
		}
	
		/*!*******************************************************************
		 * \copydoc plan::execute ()
		 *********************************************************************/
		virtual void execute () = 0;
	
	protected:
		using explicit_plan <datatype>::grid;
		using explicit_plan <datatype>::data_out;
		datatype *matrix; //!< A datatype pointer to the input data
	};
} /* one_d */

#endif /* end of include guard: PLAN_HPP_8B98W8AH */
