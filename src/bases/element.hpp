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
#include <cassert>
#include <memory>
#include "plan.hpp"
#include "solver.hpp"
#include "../utils/io.hpp"
#include "messenger.hpp"
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
	template <class datatype>
	class element
	{
	public:
		friend class plan <datatype>;
		/*!*******************************************************************
		 * \brief An iterator for the element class
		 * 
		 * This iterator steps through the scalar fields in the order they
		 * were added.
		 *********************************************************************/
		typedef std::vector <int>::iterator iterator;
		
		/*!*******************************************************************
		* \param i_name The string representation of the element
		* \param n_boundaries The integer number of boundaries (must be a multiple of 2)
		* \param i_inputParams The parameter object that contains the input parameters of the run
		* \param i_messenger_ptr A pointer to a messenger object
		* \param i_flags An integer set of execution flags
		*********************************************************************/
		element (int i_name, int n_boundaries, io::parameter_map& i_inputParams, messenger <datatype>* i_messenger_ptr, int i_flags) : inputParams (i_inputParams) {
			name = i_name;
			boundary_weights.resize (n_boundaries);
			inputParams = i_inputParams;
			messenger_ptr = i_messenger_ptr;
			flags = i_flags;
			timestep = 0.0;
			duration = 0.0;
			for (int i = 0; i < n_boundaries; ++i) {
				if (messenger_ptr->linked (i)) {
					boundary_weights [i] = 0.5;
				} else {
					boundary_weights [i] = 0.0;
				}
			}
		}
		
		virtual ~element () {}
		
		/*!*******************************************************************
		 * \brief Get the datatype reference to the named scalar
		 * 
		 * \param name The integer name from the index enumeration
		 * 
		 * This must be implemented in a subclass and will depend on the 
		 * storage system.
		 * 
		 * \return A datatype reference to the first element of the named scalar
		 *********************************************************************/
		virtual datatype& operator[] (int name) = 0;
	
		/*!*******************************************************************
		 * \brief Get the datatype reference to the given index of the named scalar
		 * 
		 * \param name The integer name from the index enumeration
		 * \param index The integer index of interest	
		 * 
		 * For simplicity, this may need to be overloaded in higher dimensions.
		 * 
		 * \return A datatype reference to the given index of the named scalar
		 *********************************************************************/
		virtual datatype& operator() (int name, int index = 0) {
			return (&((*this) [name])) [index];
		}
		
		/*!*******************************************************************
		 * \brief Get the iterator for the class at the first scalar index
		 * 
		 * This iterator has all the properties of a vector iterator. 
		 * Dereferencing it returns the scalar indices in the order they were
		 * added.
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
		 * \brief Set the collocation grid.
		 * 
		 * \param i_grid A shared_ptr to a collocation_grid object
		 * 
		 * TODO This assumes 1D n^3. Either the grid should be moved to a subclass or made more general
		 *********************************************************************/
		inline void set_grid (std::shared_ptr <collocation_grid <datatype> > i_grid) {
			grid = i_grid;
		}

		/*!*******************************************************************
		 * \brief Set the matrix solver.
		 * 
		 * \param i_solver A shared_ptr to a solver object
		 * 
		 * TODO This assumes 1 equation. It should be generalized for multiple equations.
		 *********************************************************************/
		inline void add_solver (std::shared_ptr <solver <datatype> > i_solver) {
			solvers.push_back (i_solver);
		}
		
		inline void add_name (int i_name) {
			names.push_back (i_name);
		}

		/*!*******************************************************************
		 * \brief Set the transform operation
		 * 
		 * \param i_plan A shared_ptr to the transform object
		 * 
		 * TODO This assumes one scalar field. It should be generalized.
		 *********************************************************************/
		inline void add_transform (std::shared_ptr<plan <datatype> > i_plan) {
			transforms.push_back (i_plan);
		}

		/*!*******************************************************************
		 * \brief Adds a plan to be executed in order
		 * 
		 * \param i_plan A shared pointer to the plan to add
		 *********************************************************************/
		inline void add_pre_plan (std::shared_ptr <plan <datatype> > i_plan) {
			TRACE ("Adding plan...");
			pre_transform_plans.push_back (std::move (i_plan));
			TRACE ("Added.");
		}
		
		/*!*******************************************************************
		 * \brief Adds a plan to be executed in order
		 * 
		 * \param i_plan A shared pointer to the plan to add
		 *********************************************************************/
		inline void add_post_plan (std::shared_ptr <plan <datatype> > i_plan) {
			TRACE ("Adding plan...");
			post_transform_plans.push_back (std::move (i_plan));
			TRACE ("Added.");
		}
		
		/*!*******************************************************************
		 * \brief Adds an implicit plan to be executed in order once at the start
		 * 
		 * \param i_plan A shared pointer to the plan to add
		 *********************************************************************/
		inline void add_implicit_plan (std::shared_ptr <plan <datatype> > i_plan) {
			TRACE ("Adding implicit plan...");
			implicit_plans.push_back (std::move (i_plan));
			TRACE ("Added.");
		}
		
		/*!*******************************************************************
		 * \brief Initialize the scalar name
		 * 
		 * \param name The integer name index to be initialized
		 * \param initial_conditions The datatype array of initial conditions
		 *********************************************************************/
		virtual void initialize (int name, datatype* initial_conditions = NULL) = 0;
		
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
				transform_inverse ();
			}
		}
		
		/*!*******************************************************************
		 * \brief Reset any matrices if necessary
		 * 
		 * This method should be overwritten if the element solves with an
		 * implicit part. It should only be called if the timestep duration
		 * has changed, for the most part.
		 *********************************************************************/
		virtual void implicit_reset () {
			TRACE ("Resetting implicits...");
		};
		
		/*!**********************************************************************
		 * \brief Transform from spectral space to physical space
		 * 
		 * In some cases, like the cosine transform, this can work in reverse.
		 * 
		 * TODO Multiple transforms and batch transforms should be possible
		 * TODO Need implementation if reverse transform is not forward transform
		 ************************************************************************/
		virtual void transform_inverse () {
			TRACE ("Transforming...");
			for (int i = 0; i < (int) transforms.size (); ++i) {
				transforms [i]->execute ();
			}
			if (flags & transformed) {
				flags &= ~transformed;
			} else {
				flags |= transformed;
			}
		}
		
		virtual void solve () {
			datatype t_timestep;
			t_timestep = calculate_timestep ();
			messenger_ptr->min (&t_timestep);
			
			for (int i = 0; i < (int) solvers.size (); ++i) {
				solvers [i]->execute ();
			}
			
			duration += timestep;
			INFO ("TOTAL TIME: " << duration);
			if (t_timestep != timestep) {
				flags &= ~unchanged_timestep;
				INFO ("Updating timestep: " << t_timestep);
			} else {
				flags |= unchanged_timestep;
			}
			timestep = t_timestep;
			flags |= transformed;
		}
		
		/*!**********************************************************************
		 * \brief Calculate the new timestep duration
		 * 
		 * This method should be overwritten in the final class. It uses the 
		 * knowledge of the user to beat numerical instabilities.
		 * 
		 * \return The datatype recommended timestep for the next timestep
		 ************************************************************************/
		virtual datatype calculate_timestep () = 0;
		
		/*!*******************************************************************
		 * \brief Execute the boundary conditions
		 * 
		 * In general, this should be overwritten in subclasses.
		 * 
		 * TODO I'm not entirely enthused about this method.
		 *********************************************************************/
		virtual void execute_boundaries () = 0;
		
		/*!**********************************************************************
		 * \brief The main function call of the class
		 * 
		 * This method tells the element to begin the main run of the simulation.
		 * It runs through all the specified plans in the appropriate order, and 
		 * updates the values as necessary. Output, if desired, is specified by 
		 * the output streams.
		 ************************************************************************/
		virtual void run ();
		
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
		int name; //!< An integer representation of the element, to be used in file output
		io::parameter_map& inputParams; //!< The map that contains the input parameters
		messenger <datatype>* messenger_ptr; //!< A pointer to the messenger object
		
		int flags; //!< An integer set of execution flags

		datatype duration; //!< The datatype total simulated time
		datatype timestep; //!< The datatype timestep length

		std::vector <int> names; //!< A vector of integer name indices of the contained scalars
		
		std::shared_ptr <collocation_grid <datatype> > grid; //!< A shared pointer to the collocation grid
		
		std::shared_ptr <io::output <datatype> > failsafe_dump; //!< An implementation to dump in case of failure
		std::shared_ptr <io::output <datatype> > normal_stream; //!< An implementation to output in normal space
		std::shared_ptr <io::output <datatype> > transform_stream; //!< An implementation to output in transform space

		std::vector <datatype> boundary_weights; //!< A datatype vector of boundary weights

		/*
			TODO It may make marginal more sense to move boundary_weights to the messenger class...
		*/

	private:
		std::vector<std::shared_ptr<plan <datatype> > > transforms; //!< A shared pointer to the forward transform
		std::vector<std::shared_ptr<solver <datatype> > > solvers; //!< A vector of shared pointers to the matrix solvers
		
		std::vector <std::shared_ptr <plan <datatype> > > pre_transform_plans; //!< A vector of shared pointers of plans to be executed
		std::vector <std::shared_ptr <plan <datatype> > > post_transform_plans; //!< A vector of shared pointers of plans to be executed
		std::vector <std::shared_ptr <plan <datatype> > > implicit_plans; //!< A vector of shared pointers of plans to be executed
	};
} /* bases */

#endif /* end of include guard: ELEMENT_HPP_IUTSU4TQ */
