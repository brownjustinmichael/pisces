/*!***********************************************************************
 * \file bases/solver.hpp
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
};

/*!**********************************************************************
 * \brief Flags specifying the time to execute a particular plan
 ************************************************************************/
enum solver_plan_flags {
	pre_plan = 0x01,
	mid_plan = 0x02,
	post_plan = 0x04,
	pre_solve_plan = 0x08
};

namespace plans
{
	/*!*******************************************************************
	 * \brief A class designed to solve a matrix equation
	 * 
	 * The solver class is a more specialized unit than the element class, but with broader use than the plan class. It is designed to keep track of an entire matrix equation, both in regards to its setup and solution. It is also used as a class of convenience. Since the solver class contains most of the information needed by the various plans employed by the code, it can be passed to them in lieu of the other
	 *********************************************************************/
	template <class datatype>
	class solver : public plan <datatype>
	{
	public:
		using plans::plan <datatype>::element_flags;
		using plans::plan <datatype>::component_flags;
	private:
		std::vector <std::string> deps;
	
		std::vector <std::shared_ptr <plan <datatype> > > pre_transform_plans; //!< A vector of shared pointers of plans to be executed before the transforms
		std::vector <std::shared_ptr <plan <datatype> > > mid_transform_plans; //!< A vector of shared pointers of plans to be executed after the vertical transform
		std::vector <std::shared_ptr <plan <datatype> > > post_transform_plans; //!< A vector of shared pointers of plans to be executed after both transforms
		std::vector <std::shared_ptr <plan <datatype> > > pre_solve_plans; //!< A vector of shared pointers of plans to be executed after both transforms
	
	public:
		/*!*******************************************************************
		 * \copydoc explicit_plan::explicit_plan ()
		 *********************************************************************/
		solver (int *i_element_flags, int *i_component_flags) : 
		plan <datatype> (i_element_flags, i_component_flags) {}
	
		virtual ~solver () {}
		
		/*!**********************************************************************
		 * \brief Get the version of the class
		 ************************************************************************/
		static versions::version& version () {
			static versions::version version ("1.0.1.0");
			return version;
		}
	
		virtual int n_dependencies () {
			return (int) deps.size ();
		}

		virtual const std::string& get_dependency (int i) {
			return deps [i];
		}

		virtual void add_dependency (std::string name) {
			deps.push_back (name);
		}
	
		/*!**********************************************************************
		 * \brief Return a pointer to the solver's matrix for the index dimension
		 * 
		 * \param index The index specifying from which dimension to grab the matrix
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
	
		/*!*******************************************************************
		 * \brief Solve the matrix equation
		 * 
		 * This method does not contain the actual implementation of the solution, which should be handled in _solve.
		 *********************************************************************/
		virtual void execute () = 0;
		
		class factory
		{
		public:
			virtual ~factory () {}
		
			virtual std::shared_ptr <plans::solver <datatype>> instance (grids::grid <datatype> **grids, datatype *i_data, datatype *i_rhs, int *i_element_flags = NULL, int *i_component_flags = NULL) const = 0;
		};
	};
} /* plans */

#endif /* end of include guard: SOLVER_HPP_BDH126SH */
