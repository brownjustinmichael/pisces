/*!**********************************************************************
 * \file solver_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-11.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_TWO_D_HPP_JW4BV4PS
#define SOLVER_TWO_D_HPP_JW4BV4PS

#include "../bases/messenger.hpp"
#include "../bases/solver.hpp"
#include "element_two_d.hpp"
#include "plan_two_d.hpp"

namespace two_d
{
	namespace fourier
	{
		namespace chebyshev
		{
			template <class datatype>
			class solver : public bases::solver <datatype>, explicit_plan <datatype>
			{
			public:
				solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, bases::messenger* i_messenger_ptr, datatype& i_timestep, datatype& i_alpha_0, datatype& i_aplha_n, datatype* i_data_in, datatype* i_explicit_rhs, datatype* i_implicit_rhs, datatype* i_data_out = NULL, int i_flags = 0x00);
				
				virtual ~solver () {}
				
				void execute (int &element_flags);
			
			private:
				void _factorize ();
		
				using explicit_plan <datatype>::n;
				using explicit_plan <datatype>::m;
				using explicit_plan <datatype>::data_in;
				using explicit_plan <datatype>::data_out;
				using bases::solver <datatype>::flags;
				using explicit_plan <datatype>::grid_n;
				using explicit_plan <datatype>::grid_m;

				bases::messenger* messenger_ptr;
		
				datatype& timestep; //!< A datatype reference to the current timestep
				datatype& alpha_0; //!< A datatype reference to the current edge_0 weight
				datatype& alpha_n; //!< A datatype reference to the current edge_n weight
				std::vector <datatype> values_0, values_n;

				datatype* positions;
				int temp_n;
				int excess_0; //!< The integer number of elements to recv from edge_0
				int excess_n; //!< The integer number of elements to recv from edge_n
				int ex_excess_0; //!< The integer number of elements to send to edge_0
				int ex_excess_n; //!< The integer number of elements to send to edge_n
				int ntop, nbot;
	
				datatype* explicit_rhs; //!< The datatype array of the right-hand-side of the matrix equation
				datatype* implicit_rhs; //!< The datatype array of the right-hand-side of the matrix equation
				datatype* default_matrix; //!< The datatype array of the non-timestep dependent matrix component
				datatype* matrix; //!< The datatype array of the matrix component to be timestep-multiplied
	
				std::vector <datatype> data_temp; //!< A datatype vector to be used in lieu of data_out for non-updating steps
				std::vector <datatype> positions_0; //!< A datatype vector of excess positions from edge_0
				std::vector <datatype> positions_n; //!< A datatype vector of excess positions from edge_n
				std::vector <datatype> factorized_matrix; //!< A datatype vector containing the factorized sum of default matrix and timestep * matrix
				std::vector <datatype> boundary_matrix; //!< A datatype vector containing the factorized sum of default matrix and timestep * matrix
				std::vector <datatype> previous_rhs;
				std::vector <int> ns;
				std::vector <int> ipiv; //!< A vector of integers needed to calculate the factorization
				std::vector <int> bipiv; //!< A vector of integers needed to calculate the factorization
				std::vector <datatype> horizontal_plus_matrix, horizontal_minus_matrix;
			};
		} /* chebyshev */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: SOLVER_TWO_D_HPP_JW4BV4PS */
