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
#include "../boundaries/boundary.hpp"

namespace plans
{
	namespace solvers
	{
		template <class datatype>
		class incompressible : public solver <datatype>
		{
		public:
			incompressible (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, mpi::messenger* i_messenger_ptr, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_0, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_n, datatype *i_rhs, datatype* i_data, datatype *i_data_x, datatype *i_data_z, int *i_element_flags, int *i_component_flags, int *i_component_x, int *i_component_z);
		
			virtual ~incompressible () {}
		
			datatype *matrix_ptr () {
				return NULL;
			}

			void factorize ();
			void execute ();
	
			class factory : public plans::solvers::solver <datatype>::factory
			{
			private:
				mpi::messenger *messenger_ptr;
				std::shared_ptr <boundaries::boundary <datatype>> boundary_0, boundary_n;
				plans::solvers::equation <datatype> &equation_x, &equation_z;

			public:
				factory (mpi::messenger *i_messenger_ptr, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_0, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_n, plans::solvers::equation <datatype> &i_equation_x, plans::solvers::equation <datatype> &i_equation_z) : messenger_ptr (i_messenger_ptr), boundary_0 (i_boundary_0), boundary_n (i_boundary_n), equation_x (i_equation_x), equation_z (i_equation_z) {
				}
			
				virtual ~factory () {}
			
				virtual std::shared_ptr <plans::solvers::solver <datatype>> instance (grids::grid <datatype> **grids, datatype *i_data, datatype *i_rhs, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
					return std::shared_ptr <plans::solvers::solver <datatype>> (new incompressible (*grids [0], *grids [1], messenger_ptr, boundary_0, boundary_n, i_rhs, i_data, equation_x.data_ptr (), equation_z.data_ptr (), i_element_flags, i_component_flags, equation_x.component_flags, equation_z.component_flags));
				}
			};
	
		private:
			int n;
			int ldn;
			int m;
			int count;
			int ntop, nbot;
			datatype *data, *data_x, *data_z, *new_pos;
			datatype ex_pos_0, ex_pos_m, exx_pos_0, exx_pos_m, exxx_pos_0, exxx_pos_m;
			int flags;
			int *component_flags_x, *component_flags_z;
			grids::grid <datatype> &grid_n;
			grids::grid <datatype> &grid_m;
			const datatype *pos_n;
			datatype *pos_m;
			datatype *sub_ptr, *diag_ptr, *sup_ptr;
			int excess_0, excess_n, id, np;
			datatype *rhs_ptr;
			std::shared_ptr <plans::plan <datatype> > transform, transform_h, x_deriv, z_deriv;

			mpi::messenger* messenger_ptr;
		
			std::shared_ptr <boundaries::boundary <datatype>> boundary_0, boundary_n;
		
			std::vector <datatype> x, bufferl, bufferr;
			std::vector <datatype> data_temp, positions, new_positions;
			std::vector <datatype> sup, sub, diag, supsup, matrix; //!< A datatype vector to be used in lieu of data_out for non-updating steps
			std::vector <int> ipiv, xipiv;
			using plans::solvers::solver <datatype>::element_flags;
		};
	} /* solvers */
} /* plans */

#endif /* end of include guard: SOLVER_TWO_D_HPP_JW4BV4PS */
