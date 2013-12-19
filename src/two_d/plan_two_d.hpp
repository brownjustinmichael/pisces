/*!**********************************************************************
 * \file plan_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-11.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef PLAN_TWO_D_HPP_XQ7AJI7K
#define PLAN_TWO_D_HPP_XQ7AJI7K

#include "../bases/grid.hpp"
#include "../bases/plan.hpp"
#include "../utils/utils.hpp"

namespace two_d
{
	/*!*******************************************************************
	 * \brief A subclass of plan, specific to explicit methods
	 * 
	 * These plans take input and produce output.
	 *********************************************************************/
	template <class datatype>
	class explicit_plan : public bases::plan <datatype>
	{
	public:
		/*!*******************************************************************
		 * \param i_n The integer number of elements in the data
		 * \param i_data_in The integer scalar index of the input
		 * \param i_data_out The integer scalar index of the output
		 * \copydoc plan::plan ()
		 *********************************************************************/
		explicit_plan (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out = NULL) :
		n (i_grid_n.n),
		m (i_grid_m.n),
		grid_n (i_grid_n),
		grid_m (i_grid_m),
		data_in (i_data_in),
		data_out (i_data_out ? i_data_out : i_data_in) {}

		virtual ~explicit_plan () {}

		/*!*******************************************************************
		 * \copydoc bases::plan::execute ()
		 *********************************************************************/
		virtual void execute (int &element_flags) = 0;

		int n; //!< An integer number of data elements (grid points) that collocation_1D will be built to handle
		int m;
		bases::grid <datatype> &grid_n, &grid_m;
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
		/*!*******************************************************************
		 * \param i_n The integer number of elements in a row of the square i_matrix
		 * \param i_grid A shared pointer to the collocation grid object
		 * \param i_matrix The datatype matrix to be updated
		 * \copydoc plan::plan ()
		 *********************************************************************/
		implicit_plan (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype *i_matrix_n, datatype *i_matrix_m, datatype* i_data_in, datatype* i_data_out) :
		explicit_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out), 
		matrix_n (i_matrix_n),
		matrix_m (i_matrix_m) {}
		
		/*
			TODO Take solver object instead of matrices
		*/

		virtual ~implicit_plan () {}

		/*!*******************************************************************
		 * \copydoc plan::execute ()
		 *********************************************************************/
		virtual void execute (int &element_flags) = 0;

		using explicit_plan <datatype>::grid_n;
		using explicit_plan <datatype>::grid_m;
		datatype *matrix_n; //!< A datatype pointer to the input data
		datatype *matrix_m; //!< A datatype pointer to the input data
	};
	
	template <class datatype>
	class add_scaled : public explicit_plan <datatype>
	{
	public:
		add_scaled (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype i_coeff, datatype *i_data_in, datatype *i_data_out) :
		explicit_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out),
		coeff (i_coeff) {}
		
		virtual ~add_scaled () {}
		
		virtual void execute (int &element_flags) {
			utils::matrix_add_scaled (m, 2 * (n / 2 + 1), coeff, data_in, data_out);
		}
	
	private:
		using explicit_plan <datatype>::n;
		using explicit_plan <datatype>::m;
		using explicit_plan <datatype>::data_in;
		using explicit_plan <datatype>::data_out;
		datatype coeff;
	};
} /* two_d */

#endif /* end of include guard: PLAN_TWO_D_HPP_XQ7AJI7K */
