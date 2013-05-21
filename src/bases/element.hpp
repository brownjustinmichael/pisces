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
#include "boundary.hpp"
#include "solver.hpp"
#include "transform.hpp"
#include "../utils/io.hpp"
#include "collocation.hpp"
#include "../config.hpp"

// #include "mpi.h"

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
		int flags; //!< An integer set of execution flags
	
		/*!*******************************************************************
		* \param i_flags An integer set of execution flags
		*********************************************************************/
		element (std::string i_name, io::parameter_map& i_inputParams, int i_flags) : inputParams (i_inputParams) {
			name = i_name;
			inputParams = i_inputParams;
			flags = i_flags;
			logger = config::make_logger ();
			
			n_boundaries = 0;
		}
		
		virtual ~element () {}
		/*!*******************************************************************
		 * \brief Calculate the matrix terms to be used in update
		 *********************************************************************/
		virtual void calculate ();
		
		virtual void send ();
		
		virtual void recv ();
	
		/*!*******************************************************************
		 * \brief Execute the boundary conditions
		 *********************************************************************/
		virtual void execute_boundaries ();
	
		/*!*******************************************************************
		 * \brief Update the element
		 *********************************************************************/
		virtual void update ();
		
		double get_dparam (std::string name) {
			return inputParams [name].asDouble;
		}
		
		int get_iparam (std::string name) {
			return inputParams [name].asInt;
		}
		
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
		
		virtual double& operator () (int name, int index = 0) {
			return (&((*this) [name])) [index];
		}
		
		/*!*******************************************************************
		 * \brief Reset every index < 0
		 *********************************************************************/
		virtual void explicit_reset () {
			if (!(flags & transformed)) {
				transform_forward->execute ();
			}
		}
		
		virtual void implicit_reset () {};
		
		
		inline void set_grid (std::shared_ptr<collocation_grid> i_grid) {
			grid = i_grid;
		}
		
		inline void set_solver (std::shared_ptr<plan> i_plan) {
			i_plan->associate (this);
			matrix_solver = i_plan;
		}
		
		inline void set_transform (std::shared_ptr<plan> i_plan) {
			transform_forward = i_plan;
			add_plan (transform_forward);
		}
	
		/*!*******************************************************************
		 * \brief Adds a plan to be executed in order
		 * 
		 * \param i_plan A shared pointer to the plan to add
		 *********************************************************************/
		inline void add_plan (std::shared_ptr <plan> i_plan) {
			TRACE (logger, "Adding plan..." << this);
			i_plan->associate (this);
			plans.push_back (i_plan);
			TRACE (logger, "Added.");
		}
		
		/*!*******************************************************************
		 * \brief Adds a boundary condition to execute in normal space
		 * 
		 * \param i_plan A shared pointer to the plan to add
		 *********************************************************************/
		inline void add_boundary (std::shared_ptr<boundary> i_boundary) {
			TRACE (logger, "Adding boundary...");
			++n_boundaries;
			i_boundary->associate (this);
			boundaries.push_back (std::move (i_boundary));
			TRACE (logger, "Added.");
		}
		
		virtual int get_boundary_index (int edge) = 0;
		
		virtual int get_boundary_increment (int edge) = 0;
		
		typedef std::vector <int>::iterator iterator;
		typedef std::vector <std::shared_ptr <plan>>::iterator iterator_plans;
		
		virtual iterator begin () {
			return names.begin ();
		}
		
		virtual iterator end () {
			return names.end ();
		}
		
		virtual iterator_plans begin_plans () {
			return plans.begin ();
		}
		
		virtual iterator_plans end_plans () {
			return plans.end ();
		}

		friend class plan;
		
		friend class boundary;
		friend class active_boundary;
	
	protected:
		std::string name;
		io::parameter_map& inputParams;
		
		int logger;

		double timestep; //!< The double timestep length

		std::vector <int> names;
		
		std::shared_ptr<collocation_grid> grid; //!< A shared pointer to the collocation grid
		
		std::shared_ptr<plan> transform_forward; //!< A shared pointer to the forward transform

		std::shared_ptr<io::output> failsafe_dump; //!< An implementation to dump in case of failure
		std::shared_ptr<io::output> normal_stream; //!< An implementation to output in normal space
		std::shared_ptr<io::output> transform_stream; //!< An implementation to output in transform space

	private:
		int n_boundaries; //!< The number of boundary conditions to execute
		
		std::shared_ptr<plan> matrix_solver; //!< A shared pointer to the matrix solver
		
		std::vector <std::shared_ptr <plan>> plans;
		std::vector<std::shared_ptr<boundary>> boundaries; //!< A vector of shared pointers to boundary conditions to be executed
	};
} /* bases */

#endif /* end of include guard: ELEMENT_HPP_IUTSU4TQ */
