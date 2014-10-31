/*!**********************************************************************
 * \file solver_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-11.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_TWO_D_HPP_JW4BV4PS
#define SOLVER_TWO_D_HPP_JW4BV4PS

#include "mpi/messenger.hpp"
#include "../equation.hpp"
#include "../solver.hpp"

namespace plans
{
	template <class datatype>
	class incompressible_corrector : public plans::solver <datatype>
	{
	public:
		incompressible_corrector (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, mpi::messenger* i_messenger_ptr, datatype *i_rhs, datatype* i_data, datatype *i_data_x, datatype *i_data_z, int *i_element_flags, int *i_component_flags, int *i_component_x, int *i_component_z);
		
		virtual ~incompressible_corrector () {}
		
		datatype *matrix_ptr () {
			return NULL;
		}

		void factorize ();
		void execute ();
	
		class factory : public plans::solver <datatype>::factory
		{
		private:
			mpi::messenger *messenger_ptr;
			plans::equation <datatype> &equation_x, &equation_z;

		public:
			factory (mpi::messenger *i_messenger_ptr, plans::equation <datatype> &i_equation_x, plans::equation <datatype> &i_equation_z) : messenger_ptr (i_messenger_ptr), equation_x (i_equation_x), equation_z (i_equation_z) {
			}
			
			virtual ~factory () {}
			
			virtual std::shared_ptr <plans::solver <datatype>> instance (plans::grid <datatype> **grids, datatype *i_data, datatype *i_rhs, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
				return std::shared_ptr <plans::solver <datatype>> (new incompressible_corrector (*grids [0], *grids [1], messenger_ptr, i_rhs, i_data, equation_x.data_ptr (), equation_z.data_ptr (), i_element_flags, i_component_flags, equation_x.component_flags, equation_z.component_flags));
			}
		};
	
	private:
		int n;
		int ldn;
		int m;
		int ntop, nbot;
		datatype *data, *data_x, *data_z;
		datatype ex_pos_0, ex_pos_m, exx_pos_0, exx_pos_m, exxx_pos_0, exxx_pos_m;
		int flags;
		int *component_flags_x, *component_flags_z;
		plans::grid <datatype> &grid_n;
		plans::grid <datatype> &grid_m;
		const datatype *pos_n;
		datatype *pos_m;
		datatype *sub_ptr, *diag_ptr, *sup_ptr;
		int excess_0, excess_n, id, np;
		datatype *rhs_ptr;
		std::shared_ptr <plans::plan <datatype> > transform, transform_h, x_deriv, z_deriv;

		mpi::messenger* messenger_ptr;
		
		std::vector <datatype> x, bufferl, bufferr;
		std::vector <datatype> data_temp, positions;
		std::vector <datatype> sup, sub, diag, supsup, matrix; //!< A datatype vector to be used in lieu of data_out for non-updating steps
		std::vector <int> ipiv, xipiv;
		using plans::solver <datatype>::element_flags;
	};
} /* plans */

#endif /* end of include guard: SOLVER_TWO_D_HPP_JW4BV4PS */
