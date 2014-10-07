/*!**********************************************************************
 * \file collocation_solver.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef COLLOCATION_SOLVER_HPP_76DA75C5
#define COLLOCATION_SOLVER_HPP_76DA75C5

#include "../boundary_two_d.hpp"
#include "../equation.hpp"
#include "linalg-block/solver.hpp"

namespace plans
{
	template <class datatype>
	class collocation_solver : public plans::solver <datatype>
	{
	private:
		int n;
		int ldn;
		int m;
		datatype *data;
		int flags;

		mpi::messenger* messenger_ptr;

		datatype& timestep; //!< A datatype reference to the current timestep
		datatype *rhs_ptr;

		const datatype* positions;
		int excess_0; //!< The integer number of elements to recv from edge_0
		int excess_n; //!< The integer number of elements to recv from edge_n

		datatype* default_matrix; //!< The datatype array of the non-timestep dependent matrix component

		std::vector <datatype> data_temp; //!< A datatype vector to be used in lieu of data_out for non-updating steps
		std::vector <datatype> factorized_matrix; //!< A datatype vector containing the factorized sum of default matrix and timestep * matrix
		std::vector <datatype> boundary_matrix; //!< A datatype vector containing the factorized sum of default matrix and timestep * matrix
		std::vector <datatype> previous_rhs;
		std::vector <int> ns;
		std::vector <int> ipiv; //!< A vector of integers needed to calculate the factorization
		std::vector <int> bipiv; //!< A vector of integers needed to calculate the factorization
		std::vector <datatype> matrix;
		
		std::shared_ptr <plans::boundary <datatype>> boundary_0, boundary_n;
		
		int inner_m;
		int ex_overlap_0;
		int overlap_0;
		int ex_overlap_n;
		int overlap_n;
		int lda;
		
		using plans::solver <datatype>::element_flags;
		using plans::solver <datatype>::component_flags;
		
	public:
		/*!**********************************************************************
		 * The collocation matrix is set up as 
		 * 
		 * 0 0 boundary row for above element       0 0
		 * 0 0 interpolating row for above element  0 0
		 * 0 0 [interpolating row for this element] 0 0
		 * 0 0 [boundary row for this element     ] 0 0
		 * 0 0 [matrix                            ] 0 0
		 * 0 0 [boundary row for this element     ] 0 0
		 * 0 0 [interpolating row for this element] 0 0
		 * 0 0 interpolating row for below element  0 0
		 * 0 0 boundary row for below element       0 0
		 ************************************************************************/
		collocation_solver (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, mpi::messenger* i_messenger_ptr, datatype& i_timestep, std::shared_ptr <plans::boundary <datatype>> i_boundary_0, std::shared_ptr <plans::boundary <datatype>> i_boundary_n, datatype *i_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags);
		collocation_solver (plans::equation <datatype> &i_solver, mpi::messenger* i_messenger_ptr, datatype& i_timestep, std::shared_ptr <plans::boundary <datatype>> i_boundary_0, std::shared_ptr <plans::boundary <datatype>> i_boundary_n);
		
		virtual ~collocation_solver () {}
		
		datatype *matrix_ptr () {
			return &matrix [0];
		}

		void factorize ();
		void execute ();
	};
} /* plans */

#endif /* end of include guard: COLLOCATION_SOLVER_HPP_76DA75C5 */
