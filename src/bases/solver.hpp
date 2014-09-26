/*!***********************************************************************
 * \file bases/solver.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_HPP_BDH126SH
#define SOLVER_HPP_BDH126SH
	
#include <memory>
#include "plan/plan.hpp"
#include "plan/grid.hpp"

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

namespace bases
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
		using bases::plan <datatype>::element_flags;
		using bases::plan <datatype>::component_flags;
	private:
		std::vector <int> deps;
		
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
		
		virtual int n_dependencies () {
			return (int) deps.size ();
		}

		virtual int& get_dependency (int i) {
			return deps [i];
		}

		virtual void add_dependency (int name) {
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
	};
	
	/*!*******************************************************************
	 * \brief A class designed to track and implement the solvers of a particular dataset
	 * 
	 * Note that the design of the element class expects that calling only the solve does not change the dataset. The solveed dataset must first be read back into the original for the solve to take effect.
	 *********************************************************************/
	template <class datatype>
	class master_solver
	{
	public:
		int *element_flags; //!< A pointer to the flags describing the global state of the element
		int *component_flags; //!< A pointer to the flags describing the state of the local variable
	protected:
		datatype *data;
	private:
		std::vector <std::shared_ptr <plan <datatype> > > pre_transform_plans; //!< A vector of shared pointers of plans to be executed before the transforms
		std::vector <std::shared_ptr <plan <datatype> > > mid_transform_plans; //!< A vector of shared pointers of plans to be executed after the vertical transform
		std::vector <std::shared_ptr <plan <datatype> > > post_transform_plans; //!< A vector of shared pointers of plans to be executed after both transforms
		std::vector <std::shared_ptr <plan <datatype> > > pre_solve_plans; //!< A vector of shared pointers of plans to be executed after both transforms
		
	public:
		/*!**********************************************************************
		 * \param i_element_flags A pointer to the flags describing the global state of the element
		 * \param i_component_flags A pointer to the flags describing the state of the local variable
		 ************************************************************************/
		master_solver (datatype *i_data, int *i_element_flags, int *i_component_flags) : element_flags (i_element_flags), component_flags (i_component_flags), data (i_data) {}
		/*
			TODO Have this take an array of grid pointers?
		*/	
		virtual ~master_solver () {}
		
		virtual int n_dependencies () = 0;

		virtual int& get_dependency (int i) = 0;

		virtual void add_dependency (int name, int flags = 0x00) = 0;
		
		/*!**********************************************************************
		 * \brief Return a pointer to the data associated with the solver
		 * 
		 * Each solver is implemented to solve a matrix equation for a particular variable. This returns a pointer to the first element of that variable's dataset.
		 ************************************************************************/
		virtual datatype *data_ptr () {
			return data;
		}
	
		/*!**********************************************************************
		 * \brief Return a pointer to the grid object for the index dimension
		 * 
		 * \param index The index specifying from which dimension to grab the grid
		 ************************************************************************/
		virtual grid <datatype> *grid_ptr (int index = 0) = 0;

		
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
		 * \param flags Binary flags describing the solve, from solve_flags
		 * 
		 * By default, this method will write the data into the solve class, perform the solve, and read the data back out into the element. This is chosen in order to allow for GPU usage in the future.
		 *********************************************************************/
		virtual void solve () {
			TRACE ("Executing...");
			if (!(*component_flags & factorized)) {
				factorize ();
			}
			_solve ();
			reset ();
		}
		
		virtual void add_solver (std::shared_ptr <solver <datatype> > i_solver, int flags = 0x00) = 0;
		
		virtual std::shared_ptr <bases::solver <datatype>> get_solver (int flags = 0x00) = 0;
		
		/*!*******************************************************************
		 * \brief Adds a plan to be executed
		 * 
		 * \param i_plan A pointer to the plan to add
		 * \param flags Binary flags to specify the time to execute the flag, from solver_plan_flags
		 *********************************************************************/
		inline void add_plan (std::shared_ptr <plan <datatype>> i_plan, int flags) {
			TRACE ("Adding plan...");
			if (flags & pre_plan) {
				pre_transform_plans.push_back (i_plan);
			}
			if (flags & mid_plan) {
				mid_transform_plans.push_back (i_plan);
			}
			if (flags & post_plan) {
				post_transform_plans.push_back (i_plan);
			}
			if (flags & pre_solve_plan) {
				pre_solve_plans.push_back (i_plan);
			}
			TRACE ("Added.");
		}
	
		/*!**********************************************************************
		 * \brief Execute the plans for the state given in flags
		 * 
		 * \param flags Binary flags specifying which flags should be executed, from solver_plan_flags
		 ************************************************************************/
		inline void execute_plans (int flags) {
			if (flags & pre_plan) {
				for (int i = 0; i < (int) pre_transform_plans.size (); ++i) {
					pre_transform_plans [i]->execute ();
				}
			}
			if (flags & mid_plan) {
				for (int i = 0; i < (int) mid_transform_plans.size (); ++i) {
					mid_transform_plans [i]->execute ();
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
	
	protected:
		/*!**********************************************************************
		 * \brief Transform the dataset
		 * 
		 * \params flags Binary flags that determine the type of solve (forward_vertical, forward_horizontal, inverse_vertical, inverse_horizontal)
		 * 
		 * This method contains the implementation of the solve, which must be overwritten in the subclasses.
		 ************************************************************************/
		virtual void _solve () = 0;
		
		virtual void _factorize () = 0;
	};
} /* bases */

#endif /* end of include guard: SOLVER_HPP_BDH126SH */
