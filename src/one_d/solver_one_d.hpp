/*!***********************************************************************
 * \file one_d/solver.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-13.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_HPP_N0BUAX6H
#define SOLVER_HPP_N0BUAX6H

#include <vector>
#include <memory>
#include "plan_one_d.hpp"
#include "../bases/solver.hpp"
#include "../utils/utils.hpp"

namespace bases
{
	class messenger;
} /* bases */

namespace one_d
{
	/*!*******************************************************************
	 * \brief A 1D implementation of a matrix solver
	 * 
	 * This matrix solver solves a matrix in each element individually and
	 * sends the results to the adjacent elements. This solver should be 
	 * iterated.
	 *********************************************************************/
	template <class datatype>
	class solver : public bases::solver <datatype>, public explicit_plan <datatype>
	{
	public:
		/*!*******************************************************************
		 * \param i_excess_0 The integer number of excess elements on the edge_0 side
		 * \param i_excess_n The integer number of excess elements on the edge_n side
		 * \param i_timestep A datatype reference to the current timestep
		 * \param i_alpha_0 A datatype reference to the edge_0 weight
		 * \param i_alpha_n A datatype reference to the edge_n weight
		 * \param i_default_matrix A datatype array containing the 0 order collocation matrix
		 * \param i_matrix The datatype matrix to be factorized
		 * \param i_name_rhs The integer representation of the matrix right-hand-side
		 * \copydoc bases::solver <datatype>::solver ()
		 *********************************************************************/
		solver (bases::messenger* i_messenger_ptr, int i_n, int i_n_iterations, datatype& i_timestep, datatype& i_alpha_0, datatype& i_alpha_n, datatype* positions, datatype *i_default_matrix, datatype *i_matrix, datatype* i_data_in, datatype* i_explicit_rhs, datatype* i_implicit_rhs, datatype* i_data_out = NULL, int i_flags = 0x00);
		
		virtual ~solver () {}
		
		/*!*******************************************************************
		 * \copydoc bases::solver <datatype>::execute ()
		 *********************************************************************/
		void execute ();
		
	protected:
		/*!*******************************************************************
		 * \copydoc bases::solver <datatype>::factorize ()
		 *********************************************************************/
		void _factorize ();
		
		using explicit_plan <datatype>::n;
		using explicit_plan <datatype>::data_in;
		using explicit_plan <datatype>::data_out;
		using bases::solver <datatype>::flags;

		bases::messenger* messenger_ptr;
		
		datatype& timestep; //!< A datatype reference to the current timestep
		datatype& alpha_0; //!< A datatype reference to the current edge_0 weight
		datatype& alpha_n; //!< A datatype reference to the current edge_n weight

		datatype* positions;
		int n_iterations;
		int excess_0; //!< The integer number of elements to recv from edge_0
		int excess_n; //!< The integer number of elements to recv from edge_n
		int expected_excess_0; //!< The integer number of elements to send to edge_0
		int expected_excess_n; //!< The integer number of elements to send to edge_n
		
		datatype* explicit_rhs; //!< The datatype array of the right-hand-side of the matrix equation
		datatype* implicit_rhs; //!< The datatype array of the right-hand-side of the matrix equation
		datatype* default_matrix; //!< The datatype array of the non-timestep dependent matrix component
		datatype* matrix; //!< The datatype array of the matrix component to be timestep-multiplied
		
		std::vector <datatype> error_0; //!< A datatype vector to be recved from edge_0
		std::vector <datatype> error_n; //!< A datatype vector to be recved from edge_n
		std::vector <datatype> out_error_0; //!< A datatype vector to be sent to edge_0
		std::vector <datatype> out_error_n; //!< A datatype vector to be sent to edge_n
		std::vector <datatype> data_temp; //!< A datatype vector to be used in lieu of data_out for non-updating steps
		std::vector <datatype> positions_0; //!< A datatype vector of excess positions from edge_0
		std::vector <datatype> positions_n; //!< A datatype vector of excess positions from edge_n
		std::vector <datatype> factorized_matrix; //!< A datatype vector containing the factorized sum of default matrix and timestep * matrix
		std::vector <datatype> previous_rhs;
		std::vector <int> ipiv; //!< A vector of integers needed to calculate the factorization
	};
} /* one_d */

#endif /* end of include guard: SOLVER_HPP_N0BUAX6H */
