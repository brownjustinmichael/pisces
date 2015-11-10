/*!**********************************************************************
 * \file pseudo.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-11.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef PSEUDO_SOLVER_TWO_D_HPP_JW4BV4PS
#define PSEUDO_SOLVER_TWO_D_HPP_JW4BV4PS

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
		class pseudo_incompressible : public solver
		{
		private:
			using plans::solvers::solver::data_in;
			using plans::solvers::solver::data_out;
			using plans::solvers::solver::element_flags;
			
			int n; //!< The horizontal extent of the data
			int ldn; //!< The horizontal extent of the data array
			int m; //!< The vertical extent of the data
			int ntop; //!< The extent of the overlap region on the top
			int nbot; //!< The extent of the overlap region on the bottom
			int excess_0; //!< The number of excess points are included on the top
			int excess_n; //!< The number of excess points are included on the bottom

			datatype gamma; //!< The ratio of the specific heats
			
			grids::variable &var_x; //!< The reference to the x_component of the velocity
			grids::variable &var_z; //!< The reference to the z_component of the velocity
			datatype *data; //!< A pointer to the pressure data
			datatype *data_x; //!< A pointer to the x component of the data
			datatype *data_z; //!< A pointer to the z component of the data
			datatype *pressure; //!< A pointer to the pressure data
			datatype *density; //!< A pointer to the density data
			datatype *x_vel; //!< A pointer to the x component of the velocity
			datatype *z_vel; //!< A pointer to the z component of the velocity

			int *component_flags_x; //!< The component flags for the x component of the data
			int *component_flags_z; //!< The component flags for the z component of the data
			grids::grid &grid_n; //!< A reference to the horizontal grid object
			grids::grid &grid_m; //!< A reference to the vertical grid object
			const datatype *pos_n; //!< A pointer to the horizontal position, for convenience and speed
			datatype *pos_m; //!< A pointer to the vertical position, for convenience and speed
			datatype *new_pos; //!< A pointer to the halfway points along the grid, for convenience and speed
			
			mpi::messenger* messenger_ptr; //!< A pointer to the mpi messenger object
			
			int id; //!< The mpi id of this element
			int np; //!< The total number of mpi processes
			// std::shared_ptr <plans::plan > transform, transform_h, x_deriv, z_deriv;
			
			std::shared_ptr <boundaries::boundary> boundary_0; //!< A shared pointer to the top boundary object
			std::shared_ptr <boundaries::boundary> boundary_n; //!< A shared pointer to the bottom boundary object
						
			std::vector <datatype> x; //!< Additional storage for the block banded solve


			std::vector <datatype> matrix, positions, new_positions, bufferl, bufferr, buffer, data_temp; //!< The matrix data for the banded solve
			std::vector <int> ipiv; //!< The positional swap information for the banded solve
			std::vector <int> xipiv; //!< The positional swap information for the block banded solve
			std::vector <datatype> grad_pressure, grad_density;
			std::vector <datatype> grad2_pressure;

			grids::variable &rhs;

			std::shared_ptr <plan> transform;

			datatype *oodx, *oodx2, *oodz, *oodz2;

		public:
			/*!**********************************************************************
			 * \param i_messenger_ptr A pointer to the mpi messenger
			 * \param i_boundary_0 A shared pointer to the top boundary object
			 * \param i_boundary_n A shared pointer to the bottom boundary object
			 * \param i_data A pointer to the initial data
			 * \param i_data_out A pointer to the output location for the data
			 * \param i_rhs A pointer to the right hand side of the equation
			 * \param i_data_x A reference to the x component of the velocity
			 * \param i_data_z A reference to the z component of the velocity
			 * @param i_pressure A pointer to the pressure data
			 * @param i_density A pointer to the density data
			 * @param i_gamma The specific heat ratio
			 ************************************************************************/
			pseudo_incompressible (mpi::messenger* i_messenger_ptr, std::shared_ptr <boundaries::boundary> i_boundary_0, std::shared_ptr <boundaries::boundary> i_boundary_n, grids::variable &i_data, grids::variable &i_data_out, grids::variable &i_rhs, grids::variable &i_data_x, grids::variable &i_data_z, datatype *i_density, datatype *i_pressure, datatype i_gamma = 5. / 3.);
			
			virtual ~pseudo_incompressible () {}

			/**
			 * @copydoc solver::get_state_in
			 */
			int get_state_in () {
				return real_real;
			}

			/**
			 * @copydoc solver::get_state
			 */
			int get_state () {
				return real_spectral;
			}
			
			/*!**********************************************************************
			 * \copydoc solver::matrix_ptr
			 ************************************************************************/
			datatype *matrix_ptr () {
				return NULL;
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
			class factory : public plans::solvers::solver::factory
			{
			private:
				mpi::messenger *messenger_ptr; //!< A pointer to the mpi messenger object for the solver to be constructed
				std::shared_ptr <boundaries::boundary> boundary_0; //!< A shared pointer to the top boundary for the solver to be constructed
				std::shared_ptr <boundaries::boundary> boundary_n; //!< A shared pointer to the bottom boundary for the solver to be constructed
				typename boundaries::boundary::factory *boundary_factory_0; //!< A shared pointer to the top boundary for the solver to be constructed
				typename boundaries::boundary::factory *boundary_factory_n; //!< A shared pointer to the bottom boundary for the solver to be constructed
				datatype gamma;
				plans::solvers::equation &equation_x; //!< A reference to the x-component equation
				plans::solvers::equation &equation_z; //!< A reference to the z-component equation
				datatype *pressure; //!< The pointer to the pressure data
				datatype *density; //!< The pointer to the density data

			public:
				/*!**********************************************************************
				 * \param i_messenger_ptr A pointer to the mpi messenger object for the solver to be constructed
				 * \param i_boundary_0 A shared pointer to the top boundary for the solver to be constructed
				 * \param i_boundary_n A shared pointer to the bottom boundary for the solver to be constructed
				 * \param i_equation_x A reference to the x-component equation
				 * \param i_equation_z A reference to the z-component equation
				 * @param i_density The pointer to the density data
				 * @param i_pressure The pointer to the pressure data
				 * @param i_gamma The ratio of specific heats
				 ************************************************************************/
				factory (mpi::messenger *i_messenger_ptr, std::shared_ptr <boundaries::boundary> i_boundary_0, std::shared_ptr <boundaries::boundary> i_boundary_n, plans::solvers::equation &i_equation_x, plans::solvers::equation &i_equation_z, datatype *i_density, datatype *i_pressure, datatype i_gamma = 5. / 3.) : 
				messenger_ptr (i_messenger_ptr), 
				boundary_0 (i_boundary_0), 
				boundary_n (i_boundary_n), 
				gamma (i_gamma), 
				equation_x (i_equation_x), 
				equation_z (i_equation_z), 
				pressure (i_pressure),
				density (i_density) {}

				/*!**********************************************************************
				 * \param i_messenger_ptr A pointer to the mpi messenger object for the solver to be constructed
				 * \param i_boundary_0 A boundary factory for the top boundary
				 * \param i_boundary_n A boundary factory for the bottom boundary
				 * \param i_equation_x A reference to the x-component equation
				 * \param i_equation_z A reference to the z-component equation
				 * @param i_density The pointer to the density data
				 * @param i_pressure The pointer to the pressure data
				 * @param i_gamma The ratio of specific heats
				 ************************************************************************/
				factory (mpi::messenger *i_messenger_ptr, typename boundaries::boundary::factory &i_boundary_0, typename boundaries::boundary::factory &i_boundary_n, plans::solvers::equation &i_equation_x, plans::solvers::equation &i_equation_z, datatype *i_density, datatype *i_pressure, datatype i_gamma = 5. / 3.) : 
				messenger_ptr (i_messenger_ptr), 
				boundary_factory_0 (&i_boundary_0), 
				boundary_factory_n (&i_boundary_n), 
				gamma (i_gamma), 
				equation_x (i_equation_x), 
				equation_z (i_equation_z), 
				pressure (i_pressure),
				density (i_density) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc solver::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::solvers::solver> instance (grids::variable &i_data_in, grids::variable &i_data_out, grids::variable &i_rhs) const {
					return std::shared_ptr <plans::solvers::solver> (new pseudo_incompressible (messenger_ptr, boundary_factory_0 ? boundary_factory_0->instance (i_data_in.get_grids (), false) : boundary_0, boundary_factory_n ? boundary_factory_n->instance (i_data_in.get_grids (), true) : boundary_n, i_data_in, i_data_out, i_rhs, equation_x.data_var (), equation_z.data_var (), density, pressure, gamma));
				}
			};
		};
	} /* solvers */
} /* plans */

#endif /* end of include guard: PSEUDO_SOLVER_TWO_D_HPP_JW4BV4PS */
