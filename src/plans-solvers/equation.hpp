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
		 * Note that the design of the element class expects that calling only the solve does not change the dataset. The solved dataset must first be read back into the original for the solve to take effect. In general, 
		 *********************************************************************/
		class equation
		{
		public:
			int *element_flags; //!< A pointer to the flags describing the global state of the element
			int *component_flags; //!< A pointer to the flags describing the state of the local variable
			mpi::messenger *messenger_ptr; //!< A pointer to the mpi::messenger object associated with the equation
		
		protected:
			grids::variable &data; //!< A pointer to the data held by the equation object
		
		private:
			std::vector <std::shared_ptr <plan > > plans; //!< A vector of shared pointers of plans to be executed before the transforms

		public:
			/*!**********************************************************************
			 * \param i_data A pointer to the data associated with the equation
			 * \param i_messenger_ptr A pointer to the mpi messenger object that may be needed for particular solvers
			 ************************************************************************/
			equation (grids::variable &i_data, mpi::messenger *i_messenger_ptr = NULL) : 
			element_flags (&(i_data.element_flags)), 
			component_flags (&(i_data.component_flags)), 
			messenger_ptr (i_messenger_ptr), 
			data (i_data) {
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
			 * \brief Get the state of the variable object upon output
			 *
			 * The equation solve function operates on a variable instance, which can have various states (e.g. real-real or real-spectral). The solve operation updates only one of these, leaving the variable and transformer classes to repopulate the other states. This returns the state that is updated during the evaluation of solve.
			 * 
			 * \return The state of the variable object upon output
			 ************************************************************************/
			virtual int get_state () = 0;
		
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
			virtual const equation *get_dependency (int i) = 0;
		
			/*!**********************************************************************
			 * \brief Add a dependency to one of the solvers in the equation
			 * 
			 * \param name The name of the dependency
			 * \param flags The solver flag indicating with which direction the dependency is associated
			 ************************************************************************/
			virtual void add_dependency (equation *name, int flags = 0x00) = 0;
		
			/*!**********************************************************************
			 * \brief Return a pointer to the data associated with the solver
			 * 
			 * Each solver is implemented to solve a matrix equation for a particular variable. This returns a pointer to the first element of that variable's dataset.
			 *
			 * \return The pointer to the associated data
			 ************************************************************************/
			virtual double *data_ptr () {
				return data.ptr ();
			}

			/**
			 * @brief Return a pointer to the rhs data
			 * 
			 * @param state The state of the rhs data to retrieve
			 * 
			 * @return The pointer to the rhs data
			 */
			virtual double *rhs_ptr  (int state = 0) = 0;

			/*!**********************************************************************
			 * \brief Return the variable object of the data associated with the solver
			 *
			 * If more advanced interactions with the data are needed, a variable object can be returned at increased computational cost.
			 * 
			 * \return A reference to the variable object associated with the solver
			 ************************************************************************/
			virtual grids::variable &data_var () {
				return data;
			}
			
			/*!**********************************************************************
			 * \brief Return a pointer to the grid object for the index dimension
			 * 
			 * \param index The index specifying from which dimension to grab the grid
			 ************************************************************************/
			virtual grids::grid *grid_ptr (int index = 0) = 0;
			
			/*!**********************************************************************
			 * \brief Return a pointer to the solver's matrix for the index dimension
			 * 
			 * \param index The index specifying from which dimension to grab the matrix
			 * 
			 * Note: these matrices are implementation dependent, so the implicit plans must be associated with particular matrix types.
			 ************************************************************************/
			virtual double *matrix_ptr (int index = 0) = 0;
		
			/*!**********************************************************************
			 * \brief Get a solver from the equation object
			 * 
			 * \param flags A set of integer flags describing the direction of the solve to get (x_solver, z_solver)
			 ************************************************************************/
			virtual std::shared_ptr <plans::solvers::solver> get_solver (int flags = 0x00) = 0;
	
			/*!**********************************************************************
			 * \brief Add a solver to the equation
			 * 
			 * \param i_solver A shared pointer to the solver object to add
			 * \param flags A set of integer flags describing the direction of the solve (x_solver, z_solver)
			 ************************************************************************/
			virtual void add_solver (std::shared_ptr <solver > i_solver, int flags = 0x00) = 0;
	
			/*!**********************************************************************
			 * \brief Add a solver to the equation
			 * 
			 * \param i_factory A solver factory to generate the solver object to add
			 * \param flags A set of integer flags describing the direction of the solve (x_solver, z_solver)
			 ************************************************************************/
			virtual void add_solver (const plans::solvers::solver::factory &i_factory, int flags = 0x00) = 0;
	
			/*!*******************************************************************
			 * \brief Adds a plan to be executed
			 * 
			 * \param i_plan A shared pointer to the plan to add
			 *********************************************************************/
			inline void add_plan (std::shared_ptr <plan> i_plan) {
				TRACE ("Adding plan...");
				if (!i_plan) {
					return;
				}
				plans.push_back (i_plan);
				TRACE ("Added.");
			}
	
			/*!**********************************************************************
			 * \brief Adds a plan to be executed
			 * 
			 * \param i_factory A reference to the factory from which to construct the plan
			 ************************************************************************/
			virtual void add_plan (const plan::factory &i_factory) = 0;

			/*!**********************************************************************
			 * \brief Adds a container of plans to be executed
			 * 
			 * \param i_container A reference to the factory container from which to construct the plans
			 ************************************************************************/
			virtual void add_plan (const plan::factory_container &i_container) = 0;
			
			/**
			 * @brief Setup any implicit components of the plans contained
			 * @details Each plan class is permitted the definition of a setup() member, which is called by this method. This is designed to happen whenever the base matrices (ie. the non-time-dependent parts) become invalid.
			 */
			virtual void setup_plans () {
				for (int i = 0; i < (int) plans.size (); ++i) {
					plans [i]->setup ();
				}
				*component_flags |= plans_setup;
			}
			
			/*!**********************************************************************
			 * \brief Execute the plans associated with the given flags
			 * 
			 * \param flags The flags of the plans to executed (e.g. pre_solve)
			 ************************************************************************/
			inline void execute_plans (int flags = 0x00) {
				for (int i = 0; i < (int) plans.size (); ++i) {
					if (!(flags & implicit_only) || ((flags & implicit_only) && plans [i]->implicit ())) {
						plans [i]->execute ();
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

			/**
			 * @brief Add the plan generated by a factory to this equation
			 * @details Generate an instance of a plan from the given factory and append it to the plans to be run by this equation. Warning: this does not generate a new equation, as would normally be expected of an addition operator.
			 * 
			 * @param i_factory A shared pointer to the factory from which the plan should be instantiated
			 * @return A reference to this equation with the new plan added
			 */
			equation &operator+ (const std::shared_ptr <plan::factory> i_factory) {
				return *this + plan::factory_container (i_factory);
			}

			/**
			 * @brief Add the contents of a factory container to this equation
			 * @details Append the contents of a factory container to the equation to be run as individual plans. Warning: this does not generate a new equation, as would normally be expected of an addition operator.
			 * 
			 * @param i_container The container whose contents should be added to the equation
			 * @return A reference to this equation with the new plans added
			 */
			equation &operator+ (const plan::factory_container &i_container) {
				add_plan (-1.0 * i_container);
				return *this;
			}

			/**
			 * @brief Add a constant to the equation
			 * @details Append a scalar plan to the equation. Warning: this does not generate a new equation, as would normally be expected of an addition operator.
			 * 
			 * @param scalar The scalar to add to the equation
			 * @return A reference to this equation with the new plans added
			 */
			equation &operator+ (const double scalar) {
				add_plan ((plan::factory_container) constant (-scalar));
				return *this;
			}

			/**
			 * @brief Add a new plan or plans to the equation
			 * @details Append the plans to the equation, making sure to multiply them by -1 such that the coefficients are generated properly. Warning: this does not generate a new equation, as would normally be expected of an addition operator.
			 * 
			 * @param i_other Anything that has an overloaded + operator (factory, factory_container, etc.)
			 * @return A reference to this equation with the new plan or plans added
			 */
			template <class other>
			equation &operator- (other i_other) {
				return *this + i_other * (-1.0);
			}

			/**
			 * @brief Add a new plan or plans to the equation
			 * @details Append the plans on the right side of the equal sign to the equation, making sure to multiply them by -1 such that the coefficients are generated properly. Warning: this does not generate a new equation, as would normally be expected of an addition operator.
			 * 
			 * @param i_other Anything that has an overloaded + operator (factory, factory_container, etc.)
			 * @return A reference to this equation with the new plan or plans added
			 */
			template <class other>
			equation &operator== (other i_other) {
				return *this - i_other;
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
