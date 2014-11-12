/*!**********************************************************************
 * \file collocation_solver.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef COLLOCATION_SOLVER_HPP_76DA75C5
#define COLLOCATION_SOLVER_HPP_76DA75C5

#include "linalg-block/solver.hpp"

#include "../boundary.hpp"
#include "../solver.hpp"

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
		
		virtual ~collocation_solver () {}
		
		datatype *matrix_ptr () {
			return &matrix [0];
		}

		void factorize ();
		void execute ();
		
		class factory : public plans::solver <datatype>::factory
		{
		private:
			mpi::messenger *messenger_ptr;
			datatype &timestep;
			std::shared_ptr <plans::boundary <datatype>> boundary_0, boundary_n;

		public:
			factory (mpi::messenger *i_messenger_ptr, datatype &i_timestep, std::shared_ptr <plans::boundary <datatype>> i_boundary_0 = NULL, std::shared_ptr <plans::boundary <datatype>> i_boundary_n = NULL) : messenger_ptr (i_messenger_ptr), timestep (i_timestep), boundary_0 (i_boundary_0), boundary_n (i_boundary_n) {}
			virtual ~factory () {}
			
			virtual std::shared_ptr <plans::solver <datatype>> instance (plans::grid <datatype> **grids, datatype *i_data, datatype *i_rhs, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
				return std::shared_ptr <plans::solver <datatype>> (new collocation_solver (*grids [0], *grids [1], messenger_ptr, timestep, boundary_0, boundary_n, i_rhs, i_data, i_element_flags, i_component_flags));
			}
		};
	};
} /* plans */

#endif /* end of include guard: COLLOCATION_SOLVER_HPP_76DA75C5 */
