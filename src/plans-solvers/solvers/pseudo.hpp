/*!**********************************************************************
 * \file incompressible.hpp
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
		class pseudo_incompressible : public solver <datatype>
		{
		private:
			using plans::solvers::solver <datatype>::element_flags;
			
			int n; //!< The horizontal extent of the data
			int ldn; //!< The horizontal extent of the data array
			int m; //!< The vertical extent of the data
			int ntop; //!< The extent of the overlap region on the top
			int nbot; //!< The extent of the overlap region on the bottom
			int excess_0; //!< The number of excess points are included on the top
			int excess_n; //!< The number of excess points are included on the bottom
			
			datatype *data; //!< A pointer to the pressure data
			datatype *data_x; //!< A pointer to the x component of the data
			datatype *data_z; //!< A pointer to the z component of the data
			
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
			
			std::vector <datatype> positions, diff, diff2; //!< A vector containing the vertical positions
			
			std::vector <datatype> x; //!< Additional storage for the block banded solve

			std::vector <datatype> matrix; //!< The matrix data for the banded solve
			std::vector <int> ipiv; //!< The positional swap information for the banded solve
			std::vector <int> xipiv; //!< The positional swap information for the block banded solve
			
		public:
			/*!**********************************************************************
			 * \copydoc solver::solver
			 * 
			 * \param i_grid_n The horizontal grid object
			 * \param i_grid_m The vertical grid object
			 * \param i_messenger_ptr A pointer to the mpi messenger
			 * \param i_boundary_0 A shared pointer to the top boundary object
			 * \param i_boundary_n A shared pointer to the bottom boundary object
			 * \param i_data A pointer to the initial data
			 * \param i_data_x A pointer to the x-component data
			 * \param i_data_z A pointer to the z-component data
			 * \param i_component_x A pointer to the x-component flags
			 * \param i_component_z A pointer to the z-component flags
			 ************************************************************************/
			pseudo_incompressible (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, mpi::messenger* i_messenger_ptr, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_0, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_n, datatype* i_data, datatype *i_data_x, datatype *i_data_z, int *i_element_flags, int *i_component_flags, int *i_component_x, int *i_component_z);
			
			virtual ~pseudo_incompressible () {}
			
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

				factory (mpi::messenger *i_messenger_ptr, typename boundaries::boundary <datatype>::factory &i_boundary_0, typename boundaries::boundary <datatype>::factory &i_boundary_n, plans::solvers::equation <datatype> &i_equation_x, plans::solvers::equation <datatype> &i_equation_z) : messenger_ptr (i_messenger_ptr), boundary_factory_0 (&i_boundary_0), boundary_factory_n (&i_boundary_n), equation_x (i_equation_x), equation_z (i_equation_z) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc solver::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::solvers::solver <datatype>> instance (grids::grid <datatype> **grids, datatype *i_data, datatype *i_rhs, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
					return std::shared_ptr <plans::solvers::solver <datatype>> (new pseudo_incompressible (*grids [0], *grids [1], messenger_ptr, boundary_factory_0 ? boundary_factory_0->instance (grids, false) : boundary_0, boundary_factory_n ? boundary_factory_n->instance (grids, true) : boundary_n, i_data, equation_x.data_ptr (), equation_z.data_ptr (), i_element_flags, i_component_flags, equation_x.component_flags, equation_z.component_flags));
				}
			};
		};
	} /* solvers */
} /* plans */

#endif /* end of include guard: PSEUDO_SOLVER_TWO_D_HPP_JW4BV4PS */