/*!***********************************************************************
 * \file solver.hpp
 * Spectral Element
 * 
 * This extension extends the capacities of the plans module by adding equations and solvers.
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_HPP_BDH126SH
#define SOLVER_HPP_BDH126SH
	
#include <memory>

#include "versions/version.hpp"
#include "plans/explicit_plan.hpp"
#include "plans/implicit_plan.hpp"
#include "plans/real_plan.hpp"
#include "plans/grids/grid.hpp"

/*!**********************************************************************
 * \namespace plans::solvers
 * 
 * \brief An extension of plans that contains the solver information, including the equation and solver classes.
 ************************************************************************/
namespace plans
{
	namespace solvers
	{
		/*!*******************************************************************
		 * \brief Execution flags used by the solver class
		 *********************************************************************/
		enum solver_flags {
			not_x_solver = 0x01,
			not_y_solver = 0x02,
			not_z_solver = 0x04,
			x_solver = not_y_solver | not_z_solver,
			y_solver = not_x_solver | not_z_solver,
			z_solver = not_x_solver | not_y_solver,
			implicit_rhs = 0x00,
			spectral_rhs = 0x00,
			explicit_rhs = 0x01,
			real_rhs = 0x02,
			factorized = 0x08,
			first_run = 0x100,
			ignore_net = 0x40000000
		};
		
		/*!**********************************************************************
		 * \brief Flags specifying the time to execute a particular plan
		 ************************************************************************/
		enum solver_plan_flags {
			pre_plan = 0x01,
			mid_plan = 0x02,
			post_plan = 0x04,
			pre_solve_plan = 0x08,
			implicit_only = 0x10
		};

		template <class datatype>
		class equation;
		
		/*!*******************************************************************
		 * \brief A class designed to solve a matrix equation
		 * 
		 * The solver class is a more specialized unit than the element class, but with broader use than the plan class. It is designed to keep track of an entire matrix equation, both in regards to its setup and solution. It is also used as a class of convenience. Since the solver class contains most of the information needed by the various plans employed by the code, it can be passed to them in lieu of the other
		 *********************************************************************/
		template <class datatype>
		class solver : public plan <datatype>
		{
		private:
			std::vector <equation <datatype>*> deps; //!< A vector containing the dependencies of this solver
	
			std::vector <std::shared_ptr <plan <datatype> > > pre_transform_plans; //!< A vector of shared pointers of plans to be executed before the transforms
			std::vector <std::shared_ptr <plan <datatype> > > mid_transform_plans; //!< A vector of shared pointers of plans to be executed after the vertical transform
			std::vector <std::shared_ptr <plan <datatype> > > post_transform_plans; //!< A vector of shared pointers of plans to be executed after both transforms
			std::vector <std::shared_ptr <plan <datatype> > > pre_solve_plans; //!< A vector of shared pointers of plans to be executed after both transforms
	
		public:
			using plans::plan <datatype>::element_flags;
			using plans::plan <datatype>::component_flags;
		
			/*!*******************************************************************
			 * \copydoc plan::plan ()
			 *********************************************************************/
			solver (int *i_element_flags, int *i_component_flags) : 
			plan <datatype> (i_element_flags, i_component_flags) {}
	
			virtual ~solver () {}

			virtual int type () {
				return 0;
			}
		
			/*!**********************************************************************
			 * \brief Get the version of the class
			 * 
			 * \return The version of the class
			 ************************************************************************/
			static versions::version& version () {
				static versions::version version ("1.1.0.0");
				return version;
			}
		
			/*!**********************************************************************
			 * \brief Get the number of dependencies
			 * 
			 * \return The integer number of dependencies
			 ************************************************************************/
			virtual int n_dependencies () {
				return (int) deps.size ();
			}
		
			/*!**********************************************************************
			 * \brief Get the ith dependency
			 * 
			 * \param i The index at which to get the dependency
			 * 
			 * \return The ith dependency
			 ************************************************************************/
			virtual const equation <datatype> *get_dependency (int i) {
				return deps [i];
			}
		
			/*!**********************************************************************
			 * \brief Add a dependency
			 * 
			 * \param name The name of the dependency to add
			 ************************************************************************/
			virtual void add_dependency (equation <datatype> *equation) {
				deps.push_back (equation);
			}
	
			/*!**********************************************************************
			 * \brief Return a pointer to the solver's matrix
			 * 
			 * Note: these matrices are implementation dependent, so the implicit plans must be associated with particular matrix types.
			 ************************************************************************/
			virtual datatype *matrix_ptr () = 0;
			/*
				TODO Restrict implicit plans to only be associated with one type of matrix solver
			*/
	
			/*!*******************************************************************
			 * \brief Factorize the matrix equation
			 * 
			 * This method does not check first whether the matrix has been factorized, according to the execution flags. It also does not contain the actual implementation of the factorization, which should be handled in _factorize.
			 *********************************************************************/
			virtual void factorize () = 0;
			
			virtual void setup () {}
			
			/*!*******************************************************************
			 * \brief Solve the matrix equation
			 * 
			 * This method does not contain the actual implementation of the solution, which should be handled in _solve.
			 *********************************************************************/
			virtual void execute () = 0;
		
			/*!**********************************************************************
			 * \brief An abstract factory class which equations will be able to use for convenience
			 * 
			 * Just like with the plan factories, these need to be overwritten in subclasses.
			 ************************************************************************/
			class factory
			{
			public:
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \brief Create an instance of the solver object
				 * 
				 * \param grids An array of grid objects that define the data grid
				 * \param i_data The datatype array of the data
				 * \param i_rhs The datatype array of the right hand side data
				 * \param i_element_flags A pointer to the integer flags associated with the element on the whole
				 * \param i_component_flags A pointer to the integer flags associated with the variable associated with the plan
				 * 
				 * This method creates a shared_ptr to a solver instance. The benefit to this inclusion is that the instance method can be called in a uniform way and hide communication of grid and matrix information from the user. If a plan would be created that would not do anything (e.g. something with a coefficient of 0.0), this will return a NULL shared pointer.
				 ************************************************************************/
				virtual std::shared_ptr <solver <datatype>> instance (grids::variable <datatype> &i_data, datatype *i_rhs, int *i_element_flags = NULL, int *i_component_flags = NULL) const = 0;
			};
		};
	} /* solvers */
} /* plans */

#endif /* end of include guard: SOLVER_HPP_BDH126SH */
