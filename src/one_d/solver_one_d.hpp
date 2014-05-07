/*!***********************************************************************
 * \file one_d/solver.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-13.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_HPP_N0BUAX6H
#define SOLVER_HPP_N0BUAX6H

#include "../bases/messenger.hpp"
#include <vector>
#include <memory>
#include "plan_one_d.hpp"
#include "../bases/solver.hpp"
#include "../utils/utils.hpp"

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
   	class solver : public bases::solver <datatype>
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
   		solver (bases::grid <datatype> &i_grid, bases::messenger* i_messenger_ptr, datatype& i_timestep, datatype& i_alpha_0, datatype& i_aplha_n, datatype* i_data, int *element_flags, int *component_flags);
	
   		virtual ~solver () {
   			// printf ("Destroying one_d solver\n");
   		}

   		datatype *matrix_ptr (int index = 0) {
   			return &matrix [0];
   		}
		
   		datatype *data_ptr () {
   			return data;
   		}
		
   		datatype *rhs_ptr (int index = implicit_rhs) {
   			if (index == implicit_rhs) {
   				return &implicit_rhs_vec [0];
   			} else if (index == explicit_rhs) {
   				return &explicit_rhs_vec [0];
   			} else {
   				return NULL;
   			}
   		}
		
   		bases::grid <datatype> *grid_ptr (int index = 0) {
   			return &grid;
   		}
		
   		virtual void reset () {
   			utils::scale (n, 0.0, rhs_ptr (implicit_rhs));
   			utils::scale (n, 0.0, rhs_ptr (explicit_rhs));
   		}
		
   		using bases::solver <datatype>::element_flags;
   		using bases::solver <datatype>::component_flags;
			
   	protected:
   		/*!*******************************************************************
   		 * \copydoc bases::solver <datatype>::factorize ()
   		 *********************************************************************/
   		void _factorize ();
   		void _solve ();
	
   		int n, ld;
   		bases::grid <datatype> &grid;
   		datatype *data;

   		bases::messenger* messenger_ptr;
	
   		datatype& timestep; //!< A datatype reference to the current timestep
   		datatype& alpha_0; //!< A datatype reference to the current edge_0 weight
   		datatype& alpha_n; //!< A datatype reference to the current edge_n weight
   		datatype value_0, value_n;

   		datatype* positions;
   		std::vector <datatype> matrix;
   		int temp_n;
   		int excess_0; //!< The integer number of elements to recv from edge_0
   		int excess_n; //!< The integer number of elements to recv from edge_n
   		int ex_excess_0; //!< The integer number of elements to send to edge_0
   		int ex_excess_n; //!< The integer number of elements to send to edge_n
   		int ntop, nbot;
	
   		std::vector <datatype> explicit_rhs_vec; //!< The datatype array of the right-hand-side of the matrix equation
   		std::vector <datatype> implicit_rhs_vec; //!< The datatype array of the right-hand-side of the matrix equation
   		std::vector <datatype> real_rhs_vec; //!< The datatype array of the right-hand-side of the matrix equation
   		datatype* default_matrix; //!< The datatype array of the non-timestep dependent matrix component
	
		std::vector <datatype> values_0, values_n;
   		std::vector <datatype> data_temp; //!< A datatype vector to be used in lieu of data_out for non-updating steps
   		std::vector <datatype> positions_0; //!< A datatype vector of excess positions from edge_0
   		std::vector <datatype> positions_n; //!< A datatype vector of excess positions from edge_n
   		std::vector <datatype> factorized_matrix; //!< A datatype vector containing the factorized sum of default matrix and timestep * matrix
   		std::vector <datatype> boundary_matrix; //!< A datatype vector containing the factorized sum of default matrix and timestep * matrix
   		std::vector <datatype> previous_rhs;
   		std::vector <int> ns;
   		std::vector <int> ipiv; //!< A vector of integers needed to calculate the factorization
   		std::vector <int> bipiv; //!< A vector of integers needed to calculate the factorization
   	};
} /* one_d */

#endif /* end of include guard: SOLVER_HPP_N0BUAX6H */
