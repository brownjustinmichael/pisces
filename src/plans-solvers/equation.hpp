/*!**********************************************************************
 * \file equation.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef EQUATION_HPP_2A04DA4B
#define EQUATION_HPP_2A04DA4B

#include <vector>
#include <memory>

#include "versions/version.hpp"
#include "solver.hpp"
#include "mpi/messenger.hpp"
#include "plans/source.hpp"

namespace plans
{
	namespace solvers
	{
		/*!*******************************************************************
		 * \brief A class designed to track and implement the solvers of a particular dataset
		 * 
		 * Note that the design of the element class expects that calling only the solve does not change the dataset. The solveed dataset must first be read back into the original for the solve to take effect.
		 *********************************************************************/
		template <class datatype>
		class equation
		{
		public:
			int *element_flags; //!< A pointer to the flags describing the global state of the element
			int *component_flags; //!< A pointer to the flags describing the state of the local variable
			mpi::messenger *messenger_ptr;
		
		protected:
			grids::variable <datatype> &data; //!< A pointer to the data held by the equation object
		
		private:
			std::vector <std::shared_ptr <plan <datatype> > > pre_transform_plans; //!< A vector of shared pointers of plans to be executed before the transforms
			std::vector <std::shared_ptr <plan <datatype> > > mid_transform_plans; //!< A vector of shared pointers of plans to be executed after the vertical transform
			std::vector <std::shared_ptr <plan <datatype> > > post_transform_plans; //!< A vector of shared pointers of plans to be executed after both transforms
			std::vector <std::shared_ptr <plan <datatype> > > pre_solve_plans; //!< A vector of shared pointers of plans to be executed after both transforms

		public:
			/*!**********************************************************************
			 * \param i_data A pointer to the data associated with the equation
			 * \param i_element_flags A pointer to the flags describing the global state of the element
			 * \param i_component_flags A pointer to the flags describing the state of the local variable
			 ************************************************************************/
			equation (grids::variable <datatype> &i_data, int *i_element_flags, int *i_component_flags, mpi::messenger *i_messenger_ptr = NULL) : element_flags (i_element_flags), component_flags (i_component_flags), messenger_ptr (i_messenger_ptr), data (i_data) {
				*component_flags &= ~factorized;
			}
			/*
				TODO Have this take an array of grid pointers?
			*/
		
			virtual ~equation () {}
		
			/*!**********************************************************************
			 * \brief Gets the version of the class
			 * 
			 * \return The version of the class
			 ************************************************************************/
			static versions::version& version () {
				static versions::version version ("1.0.1.0");
				return version;
			}
		
			/*!**********************************************************************
			 * \brief Get the number of dependencies in its current state
			 * 
			 * \return The number of dependencies of the equation
			 ************************************************************************/
			virtual int n_dependencies () = 0;
		
			/*!**********************************************************************
			 * \brief Get the ith dependency of the equation in its current state
			 * 
			 * \param i The index of the desired dependency
			 * 
			 * \return The ith dependency of the equation
			 ************************************************************************/
			virtual const equation <datatype> *get_dependency (int i) = 0;
		
			/*!**********************************************************************
			 * \brief Add a dependency to one of the solvers in the equation
			 * 
			 * \param name The name of the dependency
			 * \param flags The solver flag indicating with which direction the dependency is associated
			 ************************************************************************/
			virtual void add_dependency (equation <datatype> *name, int flags = 0x00) = 0;
		
			/*!**********************************************************************
			 * \brief Return a pointer to the data associated with the solver
			 * 
			 * Each solver is implemented to solve a matrix equation for a particular variable. This returns a pointer to the first element of that variable's dataset.
			 ************************************************************************/
			virtual datatype *data_ptr () {
				return data.ptr ();
			}

			virtual grids::variable <datatype> &data_var () {
				return data;
			}
			
			/*!**********************************************************************
			 * \brief Return a pointer to the grid object for the index dimension
			 * 
			 * \param index The index specifying from which dimension to grab the grid
			 ************************************************************************/
			virtual grids::grid <datatype> *grid_ptr (int index = 0) = 0;
			
			/*!**********************************************************************
			 * \brief Return a pointer to the right hand side of the matrix equation
			 * 
			 * \param index The index specifying which right hand side, (explicit_rhs, implicit_rhs, real_rhs)
			 * 
			 * The matrix solve can have several right hand sides, in particular explicit, implicit, and real. Explicit right hand sides generally contain nonlinear terms. Implicit right hand sides contain an explicit version of the implicit terms on the left hand side for additional stability. Both of these are in spectral-cartesian space for a collocation method. The real right hand side contains the explicit terms in pure Cartesian space that are transformed by the solver before a solution is calculated.
			 ************************************************************************/
			virtual datatype *rhs_ptr (int index = 0) = 0;

			/*!**********************************************************************
			 * \brief Return a pointer to the solver's matrix for the index dimension
			 * 
			 * \param index The index specifying from which dimension to grab the matrix
			 * 
			 * Note: these matrices are implementation dependent, so the implicit plans must be associated with particular matrix types.
			 ************************************************************************/
			virtual datatype *matrix_ptr (int index = 0) = 0;
		
			/*!**********************************************************************
			 * \brief Get a solver from the equation object
			 * 
			 * \param flags A set of integer flags describing the direction of the solve to get (x_solver, z_solver)
			 ************************************************************************/
			virtual std::shared_ptr <plans::solvers::solver <datatype>> get_solver (int flags = 0x00) = 0;
	
			/*!**********************************************************************
			 * \brief Add a solver to the equation
			 * 
			 * \param i_solver A shared pointer to the solver object to add
			 * \param flags A set of integer flags describing the direction of the solve (x_solver, z_solver)
			 ************************************************************************/
			virtual void add_solver (std::shared_ptr <solver <datatype> > i_solver, int flags = 0x00) = 0;
	
			/*!**********************************************************************
			 * \brief Add a solver to the equation
			 * 
			 * \param i_factory A solver factory to generate the solver object to add
			 * \param flags A set of integer flags describing the direction of the solve (x_solver, z_solver)
			 ************************************************************************/
			virtual void add_solver (const typename plans::solvers::solver <datatype>::factory &i_factory, int flags = 0x00) = 0;
	
			/*!*******************************************************************
			 * \brief Adds a plan to be executed
			 * 
			 * \param i_plan A shared pointer to the plan to add
			 * \param flags Binary flags to specify the time to execute the flag, from solver_plan_flags
			 *********************************************************************/
			inline void add_plan (std::shared_ptr <plan <datatype>> i_plan) {
				TRACE ("Adding plan...");
				if (!i_plan) {
					return;
				}
				int flags = i_plan->type ();
				if (flags & plan <datatype>::pre) {
					pre_transform_plans.push_back (i_plan);
				}
				if (flags & plan <datatype>::mid) {
					mid_transform_plans.push_back (i_plan);
				}
				if (flags & plan <datatype>::post) {
					post_transform_plans.push_back (i_plan);
				}
				if (flags & plan <datatype>::pre_solve) {
					pre_solve_plans.push_back (i_plan);
				}
				TRACE ("Added.");
			}
	
			/*!**********************************************************************
			 * \brief Adds a plan to be executed
			 * 
			 * \param i_factory A reference to the factory from which to construct the plan
			 * \param flags Binary flags to specify the time to execute the flag, from solver_plan_flags
			 ************************************************************************/
			virtual void add_plan (const typename plan <datatype>::factory &i_factory) = 0;

			virtual void add_plan (const typename plan <datatype>::factory_container &i_container) = 0;
			
			equation <datatype> &operator+ (const typename plan <datatype>::factory_container &i_container) {
				add_plan (-1.0 * i_container);
				return *this;
			}

			equation <datatype> &operator+ (const std::shared_ptr <typename plan <datatype>::factory> i_factory) {
				return *this + typename plan <datatype>::factory_container (i_factory);
			}

			equation <datatype> &operator+ (const datatype scalar) {
				add_plan (constant (-scalar));
				return *this;
			}

			template <class other>
			equation <datatype> &operator== (other i_other) {
				return *this - i_other;
			}

			template <class other>
			equation <datatype> &operator- (other i_other) {
				return *this + i_other * (-1.0);
			}

			virtual void setup_plans () {
				for (int i = 0; i < (int) pre_transform_plans.size (); ++i) {
					pre_transform_plans [i]->setup ();
				}
				for (int i = 0; i < (int) mid_transform_plans.size (); ++i) {
					mid_transform_plans [i]->setup ();
				}
				for (int i = 0; i < (int) post_transform_plans.size (); ++i) {
					post_transform_plans [i]->setup ();
				}
				for (int i = 0; i < (int) pre_solve_plans.size (); ++i) {
					pre_solve_plans [i]->setup ();
				}
				*component_flags |= plans_setup;
			}
			
			/*!**********************************************************************
			 * \brief Execute the plans associated with the given flags
			 * 
			 * \param flags The flags of the plans to executed (e.g. pre_solve)
			 ************************************************************************/
			inline void execute_plans (int flags) {
				if (flags & pre_plan) {
					for (int i = 0; i < (int) pre_transform_plans.size (); ++i) {
						pre_transform_plans [i]->execute ();
					}
				}
				if (flags & mid_plan) {
					for (int i = 0; i < (int) mid_transform_plans.size (); ++i) {
						if (!(flags & implicit_only) || ((flags & implicit_only) && mid_transform_plans [i]->implicit ())) {
							mid_transform_plans [i]->execute ();
						}
					}
				}
				if (flags & post_plan) {
					for (int i = 0; i < (int) post_transform_plans.size (); ++i) {
						post_transform_plans [i]->execute ();
					}
				}
				if (flags & pre_solve_plan) {
					for (int i = 0; i < (int) pre_solve_plans.size (); ++i) {
						pre_solve_plans [i]->execute ();
					}
				}
			}
		
			/*!**********************************************************************
			 * \brief Reset anything in the solver class that needs to be reset after the equation is solved (e.g. zero arrays, switch directions)
			 * 
			 * This method should be overwritten in subclasses if necessary.
			 ************************************************************************/
			virtual void reset () {}

			/*!*******************************************************************
			 * \brief Factorize the matrix equation
			 * 
			 * This method does not check first whether the matrix has been factorized, according to the execution flags. It also does not contain the actual implementation of the factorization, which should be handled in _factorize.
			 *********************************************************************/
			void factorize () {
				TRACE ("Factorizing...");
				_factorize ();
				*component_flags |= factorized;
			}

			/*!*******************************************************************
			 * \brief Transform the dataset according to the given flags
			 * 
			 * By default, this method will write the data into the solve class, perform the solve, and read the data back out into the element. This is chosen in order to allow for GPU usage in the future.
			 *********************************************************************/
			virtual void solve () {
				TRACE ("Executing...");
				if (!(*component_flags & plans_setup)) {
					setup_plans ();
					factorize ();
				}
				if (!(*component_flags & factorized)) {
					factorize ();
				}
				_solve ();
			}

		protected:
			/*!**********************************************************************
			 * \brief Factorize the equation
			 * 
			 * This method contains the implementation of the factorization, which must be overwritten in the subclasses.
			 ************************************************************************/
			virtual void _factorize () = 0;
		
			/*!**********************************************************************
			 * \brief Transform the dataset
			 * 
			 * This method contains the implementation of the solve, which must be overwritten in the subclasses.
			 ************************************************************************/
			virtual void _solve () = 0;
		};
	} /* solvers */
} /* plans */

#endif /* end of include guard: EQUATION_HPP_2A04DA4B */
