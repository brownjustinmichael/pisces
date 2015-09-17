/*!**********************************************************************
 * \file incompressible.hpp
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
		/*!**********************************************************************
		 * \brief A solver to correct the velocities and pressure to be incompressible
		 ************************************************************************/
		template <class datatype>
		class incompressible : public solver <datatype>
		{
		private:
			using plans::solvers::solver <datatype>::data_in;
			using plans::solvers::solver <datatype>::data_out;
			using plans::solvers::solver <datatype>::element_flags;
			
			int n; //!< The horizontal extent of the data
			int ldn; //!< The horizontal extent of the data array
			int m; //!< The vertical extent of the data
			int ntop; //!< The extent of the overlap region on the top
			int nbot; //!< The extent of the overlap region on the bottom
			int excess_0; //!< The number of excess points are included on the top
			int excess_n; //!< The number of excess points are included on the bottom
			
			grids::variable <datatype> &var_x; //!< A reference to the x component of the velocity
			grids::variable <datatype> &var_z; //!< A reference to the z component of the velocity
			datatype *data_x; //!< A pointer to the x component of the data
			datatype *data_z; //!< A pointer to the z component of the data
			datatype *new_pos; //!< A pointer to the grid positions in the solve (midpoints)
			
			datatype ex_pos_0; //!< The position one out from the top boundary
			datatype ex_pos_m; //!< The position one out from the bottom boundary
			datatype exx_pos_0; //!< The position two out from the top boundary
			datatype exx_pos_m; //!< The position two out from the bottom boundary
			datatype exxx_pos_0; //!< The position three out from the top boundary
			datatype exxx_pos_m; //!< The position three out from the bottom boundary
			
			int *component_flags_x; //!< The component flags for the x component of the data
			int *component_flags_z; //!< The component flags for the z component of the data
			grids::grid <datatype> &grid_n; //!< A reference to the horizontal grid object
			grids::grid <datatype> &grid_m; //!< A reference to the vertical grid object
			const datatype *pos_n; //!< A pointer to the horizontal position, for convenience and speed
			datatype *pos_m; //!< A pointer to the vertical position, for convenience and speed
			
			mpi::messenger* messenger_ptr; //!< A pointer to the mpi messenger object
			
			int id; //!< The mpi id of this element
			int np; //!< The total number of mpi processes
			// std::shared_ptr <plans::plan <datatype> > transform, transform_h, x_deriv, z_deriv;
			
			std::shared_ptr <boundaries::boundary <datatype>> boundary_0; //!< A shared pointer to the top boundary object
			std::shared_ptr <boundaries::boundary <datatype>> boundary_n; //!< A shared pointer to the bottom boundary object
			
			std::vector <datatype> data_temp, diff, diff2; //!< A vector that represents the actual right hand side during the solve
			std::vector <datatype> positions; //!< A vector containing the vertical positions
			std::vector <datatype> new_positions; //!< A vector containing the midpoint vertical positions
			
			std::vector <datatype> x; //!< Additional storage for the block banded solve
			std::vector <datatype> bufferl; //!< Additional storage for the block banded solve
			std::vector <datatype> bufferr; //!< Additional storage for the block banded solve
			
			std::vector <datatype> matrix; //!< The matrix data for the banded solve
			std::vector <int> ipiv; //!< The positional swap information for the banded solve
			std::vector <int> xipiv; //!< The positional swap information for the block banded solve
			
		public:
			/*!**********************************************************************
			 * \param i_messenger_ptr A pointer to the mpi messenger
			 * \param i_boundary_0 A shared pointer to the top boundary object
			 * \param i_boundary_n A shared pointer to the bottom boundary object
			 * \param i_data A pointer to the initial data
			 * \param i_data_out A pointer to the location to output the data
			 * \param i_data_x A reference to the x component of the velocity
			 * \param i_data_z A reference to the z component of the velocity
			 ************************************************************************/
			incompressible (mpi::messenger* i_messenger_ptr, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_0, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_n, grids::variable <datatype> &i_data, grids::variable <datatype> &i_data_out, grids::variable <datatype> &i_data_x, grids::variable <datatype> &i_data_z);
			
			virtual ~incompressible () {}
			
			/*!**********************************************************************
			 * \copydoc solver::matrix_ptr
			 ************************************************************************/
			datatype *matrix_ptr () {
				return NULL;
			}
		
			/**
			 * @copydoc solver::get_state
			 */
			int get_state () {
				return real_spectral;
			}
			
			/*!**********************************************************************
			 * \copydoc solver::factorize
			 ************************************************************************/
			void factorize ();
			
			/*!**********************************************************************
			 * \copydoc solver::execute
			 ************************************************************************/
			void execute ();
			
			/*!**********************************************************************
			 * \copydoc solver::factory
			 ************************************************************************/
			class factory : public plans::solvers::solver <datatype>::factory
			{
			private:
				mpi::messenger *messenger_ptr; //!< A pointer to the mpi messenger object for the solver to be constructed
				std::shared_ptr <boundaries::boundary <datatype>> boundary_0; //!< A shared pointer to the top boundary for the solver to be constructed
				std::shared_ptr <boundaries::boundary <datatype>> boundary_n; //!< A shared pointer to the bottom boundary for the solver to be constructed
				typename boundaries::boundary <datatype>::factory *boundary_factory_0; //!< A shared pointer to the top boundary for the solver to be constructed
				typename boundaries::boundary <datatype>::factory *boundary_factory_n; //!< A shared pointer to the bottom boundary for the solver to be constructed
				plans::solvers::equation <datatype> &equation_x; //!< A reference to the x-component equation
				plans::solvers::equation <datatype> &equation_z; //!< A reference to the z-component equation

			public:
				/*!**********************************************************************
				 * \param i_messenger_ptr A pointer to the mpi messenger object for the solver to be constructed
				 * \param i_boundary_0 A shared pointer to the top boundary for the solver to be constructed
				 * \param i_boundary_n A shared pointer to the bottom boundary for the solver to be constructed
				 * \param i_equation_x A reference to the x-component equation
				 * \param i_equation_z A reference to the z-component equation
				 ************************************************************************/
				factory (mpi::messenger *i_messenger_ptr, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_0, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_n, plans::solvers::equation <datatype> &i_equation_x, plans::solvers::equation <datatype> &i_equation_z) : messenger_ptr (i_messenger_ptr), boundary_0 (i_boundary_0), boundary_n (i_boundary_n), equation_x (i_equation_x), equation_z (i_equation_z) {}

				/*!**********************************************************************
				 * \param i_messenger_ptr A pointer to the mpi messenger object for the solver to be constructed
				 * \param i_boundary_0 A boundary factory for the top boundary
				 * \param i_boundary_n A boundary factory for the bottom boundary
				 * \param i_equation_x A reference to the x-component equation
				 * \param i_equation_z A reference to the z-component equation
				 ************************************************************************/
 				factory (mpi::messenger *i_messenger_ptr, typename boundaries::boundary <datatype>::factory &i_boundary_0, typename boundaries::boundary <datatype>::factory &i_boundary_n, plans::solvers::equation <datatype> &i_equation_x, plans::solvers::equation <datatype> &i_equation_z) : messenger_ptr (i_messenger_ptr), boundary_factory_0 (&i_boundary_0), boundary_factory_n (&i_boundary_n), equation_x (i_equation_x), equation_z (i_equation_z) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc solver::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::solvers::solver <datatype>> instance (grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, grids::variable <datatype> &i_rhs) const {
					return std::shared_ptr <plans::solvers::solver <datatype>> (new incompressible (messenger_ptr, boundary_factory_0 ? boundary_factory_0->instance (i_data_in.get_grids (), false) : boundary_0, boundary_factory_n ? boundary_factory_n->instance (i_data_in.get_grids (), true) : boundary_n, i_data_in, i_data_out, equation_x.data_var (), equation_z.data_var ()));
				}
			};
		};
	} /* solvers */
} /* plans */

#endif /* end of include guard: SOLVER_TWO_D_HPP_JW4BV4PS */
