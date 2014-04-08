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
#include "plan.hpp"
#include "grid.hpp"

/*!*******************************************************************
 * \brief Execution flags used by the solver class
 *********************************************************************/
enum solver_flags {
	implicit_rhs = 0x00,
	explicit_rhs = 0x01,
	real_rhs = 0x02,
	factorized = 0x08,
	first_run = 0x100,
};

namespace bases
{
	/*!*******************************************************************
	 * \brief A class designed to solve a matrix equation
	 *********************************************************************/
	template <class datatype>
	class solver : public plan <datatype>
	{
	public:
		/*!*******************************************************************
		 * \copydoc explicit_plan::explicit_plan ()
		 *********************************************************************/
		solver (int *i_element_flags, int *i_component_flags) : 
		plan <datatype> (i_element_flags, i_component_flags) {}
		
		virtual ~solver () {
			// printf ("Destroying bases solver\n");
		}
		
		/*!*******************************************************************
		 * \brief Adds a plan to be executed before either transform in order
		 * 
		 * \param i_plan A pointer to the plan to add
		 *********************************************************************/
		inline void add_pre_plan (plan <datatype>* i_plan) {
			TRACE ("Adding plan...");
			pre_transform_plans.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			TRACE ("Added.");
		}
	
		/*!*******************************************************************
		 * \brief Adds a plan to be executed after the vertical transform
		 * 
		 * \param i_plan A pointer to the plan to add
		 *********************************************************************/
		inline void add_mid_plan (plan <datatype>* i_plan) {
			TRACE ("Adding plan...");
			mid_transform_plans.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			TRACE ("Added.");
		}

		/*!*******************************************************************
		 * \brief Adds a plan to be executed after both transforms in order
		 * 
		 * \param i_plan A pointer to the plan to add
		 *********************************************************************/
		inline void add_post_plan (plan <datatype>* i_plan) {
			TRACE ("Adding plan...");
			post_transform_plans.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			TRACE ("Added.");
		}
		
		inline void add_pre_solve_plan (std::shared_ptr <plan <datatype>> i_plan) {
			TRACE ("Adding plan...");
			pre_solve_plans.push_back (i_plan);
			TRACE ("Added.");
		}

		inline void add_pre_solve_plan (plan <datatype>* i_plan) {
			TRACE ("Adding plan...");
			pre_solve_plans.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			TRACE ("Added.");
		}
		
		inline void execute_pre_plans () {
			for (int i = 0; i < (int) pre_transform_plans.size (); ++i) {
				pre_transform_plans [i]->execute ();
			}
		}
		
		inline void execute_mid_plans () {
			for (int i = 0; i < (int) mid_transform_plans.size (); ++i) {
				mid_transform_plans [i]->execute ();
			}
		}
		
		inline void execute_post_plans () {
			for (int i = 0; i < (int) post_transform_plans.size (); ++i) {
				post_transform_plans [i]->execute ();
			}
		}
		
		inline void execute_pre_solve_plans () {
			for (int i = 0; i < (int) pre_solve_plans.size (); ++i) {
				pre_solve_plans [i]->execute ();
			}
		}
		
		/*!*******************************************************************
		 * \brief Factorize the matrix equation
		 * 
		 * This method does not check first whether the matrix has been factorized, 
		 * according to the execution flags.
		 *********************************************************************/
		void factorize () {
			TRACE ("Factorizing...");
			_factorize ();
			*component_flags |= factorized;
		}
		
		virtual void reset () = 0;
		
		virtual datatype *matrix_ptr (int index = 0) = 0;
		
		virtual datatype *data_ptr () = 0;
		
		virtual datatype *rhs_ptr (int index = 0) = 0;
		
		virtual grid <datatype> *grid_ptr (int index = 0) = 0;
		
		/*!*******************************************************************
		 * \brief Solve the matrix equation
		 *********************************************************************/
		void execute () {
			TRACE ("Executing...");
			if (!(*component_flags & factorized)) {
				factorize ();
			} else {
			}
			_solve ();
			reset ();
		}
		
		using bases::plan <datatype>::element_flags;
		using bases::plan <datatype>::component_flags;
	protected:
		virtual void _factorize () = 0;
		virtual void _solve () = 0;

		std::vector <std::shared_ptr <plan <datatype> > > pre_transform_plans; //!< A vector of shared pointers of plans to be executed before the transforms
		std::vector <std::shared_ptr <plan <datatype> > > mid_transform_plans; //!< A vector of shared pointers of plans to be executed after the vertical transform
		std::vector <std::shared_ptr <plan <datatype> > > post_transform_plans; //!< A vector of shared pointers of plans to be executed after both transforms
		std::vector <std::shared_ptr <plan <datatype> > > pre_solve_plans; //!< A vector of shared pointers of plans to be executed after both transforms
	};
} /* bases */

#endif /* end of include guard: SOLVER_HPP_BDH126SH */
