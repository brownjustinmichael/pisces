/*!***********************************************************************
 * \file bases/element.hpp
 * Spectral Element
 * 
 * This file contains the element class, the basic functional unit of the
 * code.
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

namespace bases
{	
	/*!*******************************************************************
	 * \brief This is the basic class of the code
	 * 
	 * A true run will contain multiple elements linked together at the 
	 * boundaries. This code is designed to work by the collocation method.
	 *********************************************************************/
	class element
	{
	public:
		friend class plan;
		/*!*******************************************************************
		 * \brief An iterator for the element class
		 * 
		 * This iterator steps through the scalar fields in the order they
		 * were added.
		 *********************************************************************/
		typedef std::vector <int>::iterator iterator;
		
		/*!*******************************************************************
		* \param i_name The string representation of the element
		* \param i_inputParams The parameter object that contains the input parameters of the run
		* \param i_flags An integer set of execution flags
		*********************************************************************/
		element (std::string i_name, io::parameter_map& i_inputParams, int i_flags) : inputParams (i_inputParams) {
			name = i_name;
			inputParams = i_inputParams;
			flags = i_flags;
			logger = config::make_logger ();
		}
		
		virtual ~element () {}
		
		/*!*******************************************************************
		 * \brief Get the named integer parameter from the input parameter map
		 * 
		 * \param name The string name of the relevant parameter
		 * 
		 * \return The integer named parameter from the input parameter map
		 *********************************************************************/
		int get_iparam (std::string name) {
			return inputParams [name].asInt;
		}
		
		/*!*******************************************************************
		 * \brief Get the named double parameter from the input parameter map
		 * 
		 * \param name The string name of the relevant parameter
		 * 
		 * \return The double named parameter from the input parameter map
		 *********************************************************************/
		double get_dparam (std::string name) {
			return inputParams [name].asDouble;
		}
		
		/*!*******************************************************************
		 * \brief Initialize the scalar name
		 * 
		 * \param name The integer name index to be initialized
		 * \param initial_conditions The double array of initial conditions
		 *********************************************************************/
		virtual void initialize (int name, double* initial_conditions = NULL) = 0;

		/*!*******************************************************************
		 * \brief Get the double reference to the named scalar
		 * 
		 * \param name The integer name from the index enumeration
		 * 
		 * This must be implemented in a subclass and will depend on the 
		 * storage system.
		 * 
		 * \return A double reference to the first element of the named scalar
		 *********************************************************************/
		virtual double& operator[] (int name) = 0;
	
		/*!*******************************************************************
		 * \brief Get the double reference to the given index of the named scalar
		 * 
		 * \param name The integer name from the index enumeration
		 * \param index The integer index of interest	
		 * 
		 * For simplicity, this may need to be overloaded in higher dimensions.
		 * 
		 * \return A double reference to the given index of the named scalar
		 *********************************************************************/
		virtual double& operator () (int name, int index = 0) {
			return (&((*this) [name])) [index];
		}
		
		/*!*******************************************************************
		 * \brief Get the iterator for the class at the first scalar index
		 * 
		 * This iterator has all the properties of a vector iterator. 
		 * Dereferencing it returns the scalar indices in the order they were
		 * added.
		 * 
		 * TODO Change this to return them in numerical order to avoid possible ordering conflicts? Not an issue as long as all elements have the same type.
		 * 
		 * \return The iterator for the element at the first scalar index
		 *********************************************************************/
		virtual iterator begin () {
			return names.begin ();
		}
		
		/*!*******************************************************************
		 * \brief Get the iterator for the class at the last scalar index
		 * 
		 * \return The iterator for the element at the last scalar index
		 *********************************************************************/
		virtual iterator end () {
			return names.end ();
		}
		
		/*!*******************************************************************
		 * \brief Get the boundary information given an edge
		 * 
		 * \param edge The integer boundary flag for identifying the boundary
		 * \param index The integer reference to the output index of the boundary
		 * \param increment The integer reference to the output increment to the next-most inner index
		 * 
		 * To be overwritten when dimension is known.
		 *********************************************************************/
		virtual void get_boundary_info (int edge, int& index, int& increment) = 0;
		
		/*!*******************************************************************
		 * \brief Reset every scalar index < 0 and converts to spectral space
		 * 
		 * At the beginning of each timestep, all scalars with index < 0 should
		 * be reset (e.g. the right hand sides of equations have index < 0). 
		 * This must be overwritten in a subclass, which should likely call
		 * the base method. The base method converts from normal space to 
		 * spectral space if necessary.
		 *********************************************************************/
		virtual void explicit_reset () {
			if (!(flags & transformed)) {
				transform_forward->execute ();
			}
		}
	
		/*!*******************************************************************
		 * \brief Reset any matrices if necessary
		 * 
		 * This method should be overwritten if the element solves with an
		 * implicit part. It should only be called if the timestep duration
		 * has changed, for the most part.
		 *********************************************************************/
		virtual void implicit_reset () {};
	
		/*!*******************************************************************
		 * \brief Set the collocation grid.
		 * 
		 * TODO This assumes 1D n^3. Either the grid should be moved to a subclass or made more general
		 *********************************************************************/
		inline void set_grid (std::shared_ptr<collocation_grid> i_grid) {
			grid = i_grid;
		}
	
		/*!*******************************************************************
		 * \brief Set the matrix solver.
		 * 
		 * TODO This assumes 1 equation. It should be generalized for multiple equations.
		 *********************************************************************/
		inline void set_solver (std::shared_ptr<plan> i_plan) {
			i_plan->associate (this);
			matrix_solver = i_plan;
		}
	
		/*!*******************************************************************
		 * \brief Set the transform operation
		 * 
		 * TODO This assumes one scalar field. It should be generalized.
		 *********************************************************************/
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
		 * \param i_boundary A shared pointer to the boundary plan to add
		 *********************************************************************/
		inline void add_boundary (std::shared_ptr<boundary> i_boundary) {
			TRACE (logger, "Adding boundary...");
			i_boundary->associate (this);
			boundaries.push_back (std::move (i_boundary));
			TRACE (logger, "Added.");
		}
		
		/*!*******************************************************************
		 * \brief Calculate the matrix terms to be used in update
		 * 
		 * In general, this should not be overwritten in subclasses.
		 *********************************************************************/
		virtual void calculate ();
		
		/*!*******************************************************************
		 * \brief Send all relevant boundary data to adjacent elements
		 * 
		 * In general, this should not be overwritten in subclasses.
		 *********************************************************************/
		virtual void send ();
		
		/*!*******************************************************************
		 * \brief Receive all relevant boundary data from adjacent elements
		 * 
		 * In general, this should not be overwritten in subclasses.
		 *********************************************************************/
		virtual void recv ();
	
		/*!*******************************************************************
		 * \brief Execute the boundary conditions
		 * 
		 * In general, this should not be overwritten in subclasses.
		 *********************************************************************/
		virtual void execute_boundaries ();
	
		/*!*******************************************************************
		 * \brief Update the element
		 * 
		 * In general, this should not be overwritten in subclasses.
		 *********************************************************************/
		virtual void update ();

		/*!*******************************************************************
		 * \brief Output all information to a dump file in the current directory
		 * 
		 * This failsafe method is designed to create a snapshot immediately
		 * before the simulation crashed.
		 *********************************************************************/
		inline virtual void failsafe () {
			failsafe_dump->to_file ();
		}
		
	protected:
		std::string name; //!< A string representation of the element, to be used in file output
		io::parameter_map& inputParams; //!< The map that contains the input parameters
		
		int flags; //!< An integer set of execution flags
		int logger; //!< An integer representation of the logger, which is interpreted by the config class

		double timestep; //!< The double timestep length

		std::vector <int> names; //!< A vector of integer name indices of the contained scalars
		
		std::shared_ptr<collocation_grid> grid; //!< A shared pointer to the collocation grid
		
		std::shared_ptr<plan> transform_forward; //!< A shared pointer to the forward transform

		std::shared_ptr<io::output> failsafe_dump; //!< An implementation to dump in case of failure
		std::shared_ptr<io::output> normal_stream; //!< An implementation to output in normal space
		std::shared_ptr<io::output> transform_stream; //!< An implementation to output in transform space

	private:
		std::shared_ptr<plan> matrix_solver; //!< A shared pointer to the matrix solver
		
		std::vector <std::shared_ptr <plan>> plans; //!< A vector of shared pointers of plans to be executed
		std::vector<std::shared_ptr<boundary>> boundaries; //!< A vector of shared pointers to boundary conditions to be executed
	};
} /* bases */

#endif /* end of include guard: ELEMENT_HPP_IUTSU4TQ */
