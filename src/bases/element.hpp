/*!***********************************************************************
 * \file bases/element.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ELEMENT_HPP_IUTSU4TQ
#define ELEMENT_HPP_IUTSU4TQ

#include <string>
#include "plan.hpp"
#include "solver.hpp"
#include "../utils/io.hpp"
#include "collocation.hpp"
#include "../config.hpp"

namespace bases
{	
/*!*******************************************************************
 * \brief This is the basic class of the code
 * 
 * A true run will contain multiple elements tied together at the 
 * boundaries.
 *********************************************************************/
class element
{
public:
		/*!*******************************************************************
		* \param i_flags An integer set of execution flags
		*********************************************************************/
		element (std::string i_name, int i_flags) {
			name = i_name;
			flags = i_flags;
			logger = config::make_logger ();
		
			n_explicit_grid_plans = 0;
			n_explicit_space_plans = 0;
			n_implicit_plans = 0;
			n_boundaries = 0;
		
			previous_timestep = 0.0;
		}
		
		virtual ~element () {}
		/*!*******************************************************************
		 * \brief Calculate the matrix terms to be used in update
		 *********************************************************************/
		virtual void calculate ();
	
		/*!*******************************************************************
		 * \brief Execute the boundary conditions
		 *********************************************************************/
		virtual void execute_boundaries ();
	
		/*!*******************************************************************
		 * \brief Update the element
		 *********************************************************************/
		virtual void update ();
		
		inline virtual void failsafe () {
			failsafe_dump->to_file ();
		}
		
		/*!*******************************************************************
		 * \brief Get the double pointer to the named scalar
		 * 
		 * \param name The integer name from the index enumeration
		 * 
		 * \return A double pointer to the first element of the named scalar
		 *********************************************************************/
		virtual double& operator[] (int name) = 0;
		
		/*!*******************************************************************
		 * \brief Reset every index < 0
		 *********************************************************************/
		virtual void reset () {
			if (previous_timestep != timestep) {
				flags &= ~factorized;
			}
		}
		
		inline void set_grid (std::shared_ptr<collocation_grid> i_grid) {
			grid = i_grid;
		}
		
		inline void set_solver (std::shared_ptr<solver> i_solver) {
			matrix_solver = i_solver;
		}
		
		inline void set_timestep (std::shared_ptr<plan> i_plan) {
			i_plan->associate (this);
			timestep_plan = i_plan;
		}
		
		inline void set_transform (std::shared_ptr<plan> i_plan) {
			i_plan->associate (this);
			transform_forward = i_plan;
		}
	
		/*!*******************************************************************
		 * \brief Adds an explicit plan to be evaluated in grid space
		 * 
		 * \param i_plan A shared pointer to the plan to add
		 *********************************************************************/
		inline void add_explicit_grid_plan (std::shared_ptr<plan> i_plan) {
			TRACE (logger, "Adding explicit grid plan...");
			++n_explicit_grid_plans;
			i_plan->associate (this);
			explicit_grid_plans.push_back (std::move (i_plan));
			TRACE (logger, "Added.");
		}

		/*!*******************************************************************
		 * \brief Adds an explicit plan to be evaluated in normal space
		 * 
		 * \param i_plan A shared pointer to the plan to add
		 *********************************************************************/	
		inline void add_explicit_space_plan (std::shared_ptr<plan> i_plan) {
			TRACE (logger, "Adding explicit space plan...");
			++n_explicit_space_plans;
			i_plan->associate (this);
			explicit_space_plans.push_back (std::move (i_plan));
			TRACE (logger, "Added.");
		}
	
		/*!*******************************************************************
		 * \brief Adds an implicit plan to be evaluated
		 * 
		 * \param i_plan A shared pointer to the plan to add
		 *********************************************************************/
		inline void add_implicit_plan (std::shared_ptr<plan> i_plan) {
			TRACE (logger, "Adding implicit plan...");
			++n_implicit_plans;
			i_plan->associate (this);
			implicit_plans.push_back (std::move (i_plan));
			TRACE (logger, "Added.");
		}
	
		/*!*******************************************************************
		 * \brief Adds a boundary condition to execute in normal space
		 * 
		 * \param i_plan A shared pointer to the plan to add
		 *********************************************************************/
		inline void add_boundary (std::shared_ptr<plan> i_plan) {
			TRACE (logger, "Adding boundary...");
			++n_boundaries;
			i_plan->associate (this);
			boundaries.push_back (std::move (i_plan));
			TRACE (logger, "Added.");
		}
		
		friend class plan;
	
	protected:
		std::string name;
		
		int flags; //!< An integer set of execution flags
		int logger;

		double timestep; //!< The double timestep length
		
		std::shared_ptr<collocation_grid> grid; //!< A shared pointer to the collocation grid
		
		std::shared_ptr<plan> transform_forward; //!< A shared pointer to the forward transform

		std::shared_ptr<io::output> failsafe_dump; //!< An implementation to dump in case of failure
		std::shared_ptr<io::output> normal_stream; //!< An implementation to output in normal space
		std::shared_ptr<io::output> transform_stream; //!< An implementation to output in transform space

	private:
		int n_explicit_grid_plans; //!< The number of explicit grid plans to execute
		int n_explicit_space_plans; //!< The number of explicit space plans to execute
		int n_implicit_plans; //!< The number of implicit plans to execute
		int n_boundaries; //!< The number of boundary conditions to execute
		
		double previous_timestep; //!< The double duration of the previous timestep
	
		std::shared_ptr<plan> timestep_plan;
		std::shared_ptr<solver> matrix_solver; //!< A shared pointer to the matrix solver
		
		std::vector<std::shared_ptr<plan>> explicit_grid_plans; //!< A vector of shared pointers to explicit grid plans to be executed
		std::vector<std::shared_ptr<plan>> explicit_space_plans; //!< A vector of shared pointers to explicit space plans to be executed
		std::vector<std::shared_ptr<plan>> implicit_plans; //!< A vector of shared pointers to implicit plans to be executed
		std::vector<std::shared_ptr<plan>> boundaries; //!< A vector of shared pointers to boundary conditions to be executed
	};
} /* bases */

#endif /* end of include guard: ELEMENT_HPP_IUTSU4TQ */
