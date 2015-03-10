/*!**********************************************************************
 * \file laplace_solver.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef LAPLACE_SOLVER_HPP_2AB7E5D4
#define LAPLACE_SOLVER_HPP_2AB7E5D4

#include "mpi/messenger.hpp"
#include "plans/plan.hpp"

#include "../solver.hpp"
#include "../boundary.hpp"

namespace plans
{
	template <class datatype>
	class laplace_solver : public plans::solver <datatype>
	{
	public:
		laplace_solver (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, mpi::messenger* i_messenger_ptr, datatype *i_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags);
		
		virtual ~laplace_solver () {}
		
		datatype *matrix_ptr () {
			return NULL;
		}

		void factorize ();
		void execute ();
	
	private:
		int n;
		int ldn;
		int m;
		datatype *data;
		datatype ex_pos_0, ex_pos_m;
		int flags;
		grids::grid <datatype> &grid_n;
		grids::grid <datatype> &grid_m;
		const datatype *pos_n, *pos_m;
		datatype *sub_ptr, *diag_ptr, *sup_ptr;
		int excess_0, excess_n, id, np;
		datatype *rhs_ptr;

		mpi::messenger* messenger_ptr;
		
		std::vector <datatype> x;
		std::vector <datatype> sup, sub, diag, supsup; //!< A datatype vector to be used in lieu of data_out for non-updating steps
		std::vector <int> ipiv, xipiv;
	};
} /* plans */

#endif /* end of include guard: LAPLACE_SOLVER_HPP_2AB7E5D4 */
