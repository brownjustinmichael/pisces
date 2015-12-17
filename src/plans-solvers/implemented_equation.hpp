/*!**********************************************************************
 * \file implemented_equation.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef IMPLEMENTED_EQUATION_HPP_74410EFA
#define IMPLEMENTED_EQUATION_HPP_74410EFA

#include "equation.hpp"
#include "plans-transforms/transform.hpp"
 #include "plans-transforms/implemented_transformer.hpp"

namespace plans
{
	namespace solvers
	{
		/*!**********************************************************************
		 * \brief An implementation of the equation class in 2D
		 * 
		 * \copydoc equation
		 * 
		 * This implements the equation class into a two-dimensional object. It should be possible to remove this class by making equation more general.
		 ************************************************************************/
		class implemented_equation : public equation
		{
		private:
			using plans::solvers::equation::data;
			using plans::solvers::equation::element_flags;
			using plans::solvers::equation::component_flags;
			
			int n; //!< The integer number of elements in the horizontal
			int ldn; //!< The integer length of the array in the horizontal (can be larger than n)
			int m; //!< The integer number of elements in the vertical
			int flags; //!< The integer flags associated with the equation
			grids::grid &grid_n; //!< The grid object in the horizontal
			grids::grid &grid_m; //!< The grid object in the vertical
			
			std::shared_ptr <plans::solvers::solver > x_solver; //!< A pointer to the horizontal solver
			std::shared_ptr <plans::solvers::solver > z_solver; //!< A pointer to the vertical solver
			
			std::shared_ptr <grids::variable> old_rhs_ptr; //!< A variable of the old rhs pointer, for accuracy
			std::shared_ptr <grids::variable> old2_rhs_ptr; //!< A variable of the rhs pointer from two steps ago
			std::shared_ptr <grids::variable> old3_rhs_ptr; //!< A variable of the rhs pointer from three steps ago
			std::shared_ptr <grids::variable> new_rhs_ptr; //!< A variable of the current rhs pointer
			std::shared_ptr <grids::variable> cor_rhs_ptr; //!< A variable of the corrected rhs pointer, taking into account the older rhs pointers
			
			std::shared_ptr <plans::plan > transform; //!< A shared pointer to a transform if the equation has a real_rhs_vec
			std::shared_ptr <plans::transforms::transformer> transformer;

		public:
			/*!**********************************************************************
			 * \copydoc equation::equation
			 ************************************************************************/
			implemented_equation (grids::variable &i_data, mpi::messenger *i_messenger_ptr) : 
			plans::solvers::equation (i_data, i_messenger_ptr), 
			n (i_data.get_grid (0).get_n ()), 
			ldn (i_data.get_grid (0).get_ld ()), 
			m (i_data.get_grid (1).get_n ()), 
			grid_n (i_data.get_grid (0)), 
			grid_m (i_data.get_grid (1)) {
				flags = 0x00;
				new_rhs_ptr = std::shared_ptr <grids::variable> (new grids::variable (grid_n, grid_m, *element_flags));
				old_rhs_ptr = std::shared_ptr <grids::variable> (new grids::variable (grid_n, grid_m, *element_flags));
				old2_rhs_ptr = std::shared_ptr <grids::variable> (new grids::variable (grid_n, grid_m, *element_flags));
				old3_rhs_ptr = std::shared_ptr <grids::variable> (new grids::variable (grid_n, grid_m, *element_flags));
				cor_rhs_ptr = std::shared_ptr <grids::variable> (new grids::variable (grid_n, grid_m, *element_flags));
				transform = std::shared_ptr <plans::plan > (new plans::transforms::horizontal (*new_rhs_ptr, real_real));
				transformer = std::shared_ptr <plans::transforms::transformer> (new plans::transforms::implemented_transformer (i_data, component_flags, element_flags));
			}
			
			virtual ~implemented_equation () {}
			
			/*!**********************************************************************
			 * \copydoc equation::n_dependencies
			 ************************************************************************/
			virtual int n_dependencies () {
				if (*component_flags & x_solve) {
					if (x_solver) {
						return x_solver->n_dependencies ();
					} else {
						return 0;
					}
				} else if (*component_flags & z_solve) {
					if (z_solver) {
						return z_solver->n_dependencies ();
					} else {
						return 0;
					}
				} else {
					return 0;
				}
			}
			
			/*!**********************************************************************
			 * \copydoc equation::get_dependency
			 ************************************************************************/
			virtual const equation *get_dependency (int i) {
				if (*component_flags & x_solve) {
					return x_solver->get_dependency (i);
				} else if (*component_flags & z_solve) {
					return z_solver->get_dependency (i);
				} else {
					FATAL ("Cannot get dependency without valid solve direction " << *component_flags);
					throw 0;
				}
			}

			/**
			 * @brief Get the state that the equation will output to
			 * @details Each solver outputs to one state in a variable. The equation class will progress through several solvers and output at last in one state in the data. This returns that state.
			 * @return The integer state that the equation outputs to
			 */
			int get_state () {
				if (z_solver) {
					DEBUG ("Z Solve state is " << z_solver->get_state ());
					return z_solver->get_state ();
				} else if (x_solver) {
					DEBUG ("X Solve state is " << x_solver->get_state ());
					return x_solver->get_state ();
				}
				return 0;
			}

			/**
			 * @brief Return a pointer to the rhs data
			 * 
			 * @param state The state of the rhs data to retrieve
			 * 
			 * @return The pointer to the rhs data
			 */
			virtual double *rhs_ptr  (int state = 0) {
				return new_rhs_ptr->ptr (state);
			}
			
			/*!**********************************************************************
			 * \copydoc equation::get_solver
			 ************************************************************************/
			virtual std::shared_ptr <plans::solvers::solver> get_solver (int flags = 0x00) {
				if (!(flags & not_x_solver)) {
					return x_solver;
				}
				if (!(flags & not_z_solver)) {
					return z_solver;
				}
				FATAL ("Invalid flags: " << flags);
				throw 0;
			}
			
			/*!**********************************************************************
			 * \copydoc equation::add_dependency
			 ************************************************************************/
			virtual void add_dependency (equation *name, int flags = 0x00) {
				if (!(flags & not_x_solver)) {
					if (x_solver) {
						x_solver->add_dependency (name);
					}
				} else if (!(flags & not_z_solver)) {
					if (z_solver) {
						z_solver->add_dependency (name);
					}
				}
			}
			
			/*!**********************************************************************
			 * \copydoc equation::grid_ptr
			 ************************************************************************/
			grids::grid *grid_ptr (int index = 0) {
				if (index == 0) {
					return &grid_n;
				} else {
					return &grid_m;
				}
			}
			
			/*!**********************************************************************
			 * \copydoc equation::matrix_ptr
			 ************************************************************************/
			double *matrix_ptr (int index = 0) {
				if (index == 0) {
					if (x_solver) {
						return x_solver->matrix_ptr ();
					} else {
						return NULL;
					}
				} else {
					if (z_solver) {
						return z_solver->matrix_ptr ();
					} else {
						return NULL;
					}
				}
			}
			
			/*!**********************************************************************
			 * \copydoc equation::reset
			 ************************************************************************/
			virtual void reset () {
				linalg::scale (ldn * m, 0.0, new_rhs_ptr->ptr (real_real));
				linalg::scale (ldn * m, 0.0, new_rhs_ptr->ptr (real_spectral));
				linalg::scale (m * ldn, 0.0, new_rhs_ptr->ptr (spectral_spectral));
			}


			/*!**********************************************************************
			 * \copydoc equation::add_solver
			 ************************************************************************/
			virtual void add_solver (std::shared_ptr <plans::solvers::solver> i_solver, int flags = 0x00) {
				TRACE ("Adding solver...");
				if (!(flags & not_x_solver)) {
					x_solver = i_solver;
				}
				if (!(flags & not_z_solver)) {
					z_solver = i_solver;
				}
				TRACE ("Added.");
			}
			
			/*!**********************************************************************
			 * \copydoc equation::add_solver(const typename solver<double>::factory&,int)
			 ************************************************************************/
			virtual void add_solver (const plans::solvers::solver::factory &i_factory, int flags = 0x00) {
				plans::solvers::implemented_equation::add_solver (i_factory.instance (data, data, *cor_rhs_ptr), flags);
			}
			
			/*!**********************************************************************
			 * \copydoc equation::add_plan(const typename plan<double>::factory&)
			 ************************************************************************/
			void add_plan (const plans::plan::factory &i_factory) {
				TRACE ("Adding plan...");
				double* matrices [2] = {matrix_ptr (0), matrix_ptr (1)};
				plans::solvers::equation::add_plan (i_factory.instance (matrices, data, *new_rhs_ptr));
			}
			
			/*!**********************************************************************
			 * \copydoc equation::add_plan(const typename plan<double>::factory_container&)
			 ************************************************************************/
			void add_plan (const plans::plan::factory_container &i_container) {
				TRACE ("Adding plan...");
				for (int i = 0; i < (int) i_container.facts.size (); ++i)
				{
					if (i_container.facts [i]) add_plan (*(i_container.facts [i]));
				}
			}

			/**
			 * @copydoc equation::setup_plans
			 */
			void setup_plans () {
				if (x_solver) x_solver->setup ();
				if (z_solver) z_solver->setup ();
				equation::setup_plans ();
			}
		
		protected:
			/*!**********************************************************************
			 * \copydoc equation::_factorize
			 ************************************************************************/
			virtual void _factorize () {
				if (x_solver) {
					x_solver->factorize ();
				}
				if (z_solver && (x_solver != z_solver)) {
					z_solver->factorize ();
				}
			}
			
			/*!**********************************************************************
			 * \copydoc equation::_solve
			 ************************************************************************/
			virtual void _solve ();
		};
	} /* solvers */
} /* plans */

#endif /* end of include guard: IMPLEMENTED_EQUATION_HPP_74410EFA */
