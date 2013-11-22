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

#include "messenger.hpp"
#include <string>
#include <cassert>
#include <memory>
#include "plan.hpp"
#include "solver.hpp"
#include "../utils/io.hpp"
#include "grid.hpp"
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
		/*!*******************************************************************
		* \param i_name The string representation of the element
		* \param n_boundaries The integer number of boundaries (must be a multiple of 2)
		* \param i_params The parameter object that contains the input parameters of the run
		* \param i_messenger_ptr A pointer to a messenger object
		* \param i_flags An integer set of execution flags
		*********************************************************************/
		element (int i_name, int i_dimensions, io::parameters <datatype>& i_params, messenger* i_messenger_ptr, int i_flags) : 
		params (i_params),
		flags (i_flags) {
			name = i_name;
			grids.resize (i_dimensions);
			messenger_ptr = i_messenger_ptr;
			timestep = 0.0;
			duration = 0.0;
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
		inline datatype& operator[] (int name) {
			if (scalars.find (name) == scalars.end ()) {
				FATAL ("Index " << name << " not found in element.");
				throw 0;
			}
			return scalars [name] [0];
		}
		
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
		
		virtual datatype* pointer (int name, int index = 0) {
			return &((*this) [name]) + index;
		}
		
		/*!*******************************************************************
		 * \brief Set the collocation grid.
		 * 
		 * \param i_grid A shared_ptr to a grid object
		 * 
		 * TODO This assumes 1D n^3. Either the grid should be moved to a subclass or made more general
		 *********************************************************************/
		inline void set_grid (grid <datatype>* i_grid, int index = 0) {
			grids [index] = std::shared_ptr <grid <datatype>> (i_grid);
		}

		/*!*******************************************************************
		 * \brief Set the matrix solver.
		 * 
		 * \param i_solver A shared_ptr to a solver object
		 * 
		 * TODO This assumes 1 equation. It should be generalized for multiple equations.
		 *********************************************************************/
		inline void add_solver (solver <datatype>* i_solver) {
			solvers.push_back (std::shared_ptr <solver <datatype>> (i_solver));
		}

		/*!*******************************************************************
		 * \brief Set the transform operation
		 * 
		 * \param i_plan A shared_ptr to the transform object
		 * 
		 * TODO This assumes one scalar field. It should be generalized.
		 *********************************************************************/
		inline void add_forward_horizontal_transform (plan <datatype>* i_plan) {
			forward_horizontal_transforms.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			// transforms.push_back (i_plan));
		}
		
		inline void add_inverse_horizontal_transform (plan <datatype>* i_plan) {
			inverse_horizontal_transforms.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			// transforms.push_back (i_plan));
		}

		inline void add_forward_vertical_transform (plan <datatype>* i_plan) {
			forward_vertical_transforms.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			// transforms.push_back (i_plan));
		}
		
		inline void add_inverse_vertical_transform (plan <datatype>* i_plan) {
			inverse_vertical_transforms.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			// transforms.push_back (i_plan));
		}

		/*!*******************************************************************
		 * \brief Adds a plan to be executed in order
		 * 
		 * \param i_plan A shared pointer to the plan to add
		 *********************************************************************/
		inline void add_pre_plan (plan <datatype>* i_plan) {
			TRACE ("Adding plan...");
			pre_transform_plans.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			TRACE ("Added.");
		}
		
		inline void add_mid_plan (plan <datatype>* i_plan) {
			TRACE ("Adding plan...");
			mid_transform_plans.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			TRACE ("Added.");
		}

		/*!*******************************************************************
		 * \brief Adds a plan to be executed in order
		 * 
		 * \param i_plan A shared pointer to the plan to add
		 *********************************************************************/
		inline void add_post_plan (plan <datatype>* i_plan) {
			TRACE ("Adding plan...");
			post_transform_plans.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			TRACE ("Added.");
		}
		
		/*!*******************************************************************
		 * \brief Initialize the scalar name
		 * 
		 * \param name The integer name index to be initialized
		 * \param initial_conditions The datatype array of initial conditions
		 *********************************************************************/
		virtual void initialize (int name, datatype* initial_conditions = NULL, int flags = 0x00) = 0;
		
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
			TRACE ("Resetting explicits...");
			if (!(flags & transformed_horizontal)) {
				transform_horizontal_forward ();
			}
			if (!(flags & transformed_vertical)) {
				transform_vertical_forward ();
			}
		}
		
		/*!**********************************************************************
		 * \brief Transform from spectral space to physical space
		 * 
		 * In some cases, like the cosine transform, this can work in reverse.
		 * 
		 * TODO Multiple transforms and batch transforms should be possible
		 * TODO Need implementation if reverse transform is not forward transform
		 ************************************************************************/
		virtual void transform_horizontal_forward () {
			TRACE ("Transforming...");
			if (!(flags & transformed_horizontal)) {
				for (int i = 0; i < (int) forward_horizontal_transforms.size (); ++i) {
					forward_horizontal_transforms [i]->execute (flags);
				}
				flags |= transformed_horizontal;
			}
		}
		
		virtual void transform_horizontal_inverse () {
			TRACE ("Inverting...");
			if (flags & transformed_horizontal) {
				for (int i = 0; i < (int) inverse_horizontal_transforms.size (); ++i) {
					inverse_horizontal_transforms [i]->execute (flags);
				}
				flags &= ~transformed_horizontal;
			}
		}
		
		virtual void transform_vertical_forward () {
			TRACE ("Transforming...");
			if (!(flags & transformed_vertical)) {
				for (int i = 0; i < (int) forward_vertical_transforms.size (); ++i) {
					forward_vertical_transforms [i]->execute (flags);
				}
				flags |= transformed_vertical;
			}
		}
		
		virtual void transform_vertical_inverse () {
			TRACE ("Inverting...");
			if (flags & transformed_vertical) {
				for (int i = 0; i < (int) inverse_vertical_transforms.size (); ++i) {
					inverse_vertical_transforms [i]->execute (flags);
				}
				flags &= ~transformed_vertical;
			}
		}
		
		virtual void output () {
			if (normal_stream) {
				normal_stream->to_file ();
			}
			if (transform_stream) {
				transform_stream->to_file ();
			}
		}

		virtual void factorize () {
			if (!(flags & factorized)) {
				for (int i = 0; i < (int) solvers.size (); ++i) {
					solvers [i]->factorize ();
				}
				flags |= factorized;
			}
		}
		
		virtual void solve () {
			TRACE ("Beginning solve...");
			datatype t_timestep;
			t_timestep = calculate_timestep ();
			messenger_ptr->min (&t_timestep);
			
			for (int i = 0; i < (int) solvers.size (); ++i) {
				solvers [i]->execute (flags);
			}
			
			duration += timestep;
			INFO ("TOTAL TIME: " << duration);
			if (t_timestep != timestep) {
				flags &= ~unchanged_timestep;
				flags &= ~factorized;
				INFO ("Updating timestep: " << t_timestep);
			} else {
				flags |= unchanged_timestep;
			}
			timestep = t_timestep;
			if (!(flags & transformed_vertical)) {
				transform_vertical_forward ();
			}
			TRACE ("Solve complete.");
		}
		
		/*!**********************************************************************
		 * \brief Calculate the new timestep duration
		 * 
		 * This method should be overwritten in the final class. It uses the 
		 * knowledge of the user to beat numerical instabilities.
		 * 
		 * \return The datatype recommended timestep for the next timestep
		 ************************************************************************/
		virtual datatype calculate_timestep () {
			return (datatype) 0.0;
		}
		
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
		io::parameters <datatype>& params; //!< The map that contains the input parameters
		messenger* messenger_ptr; //!< A pointer to the messenger object
		
		int flags; //!< An integer set of execution flags

		datatype duration; //!< The datatype total simulated time
		datatype timestep; //!< The datatype timestep length

		std::map <int, std::vector <datatype> > scalars; //!< A vector of scalar vectors
		std::vector <std::shared_ptr <grid <datatype> > > grids; //!< A shared pointer to the collocation grid
		
		std::shared_ptr <io::output> failsafe_dump; //!< An implementation to dump in case of failure
		std::shared_ptr <io::output> normal_stream; //!< An implementation to output in normal space
		std::shared_ptr <io::output> transform_stream; //!< An implementation to output in transform space

	private:
		// std::vector<plan <datatype>* > transforms; //!< A shared pointer to the forward transform
		std::vector<std::shared_ptr<plan <datatype> > > forward_horizontal_transforms; //!< A shared pointer to the forward transform
		std::vector<std::shared_ptr<plan <datatype> > > inverse_horizontal_transforms; //!< A shared pointer to the forward transform
		std::vector<std::shared_ptr<plan <datatype> > > forward_vertical_transforms; //!< A shared pointer to the forward transform
		std::vector<std::shared_ptr<plan <datatype> > > inverse_vertical_transforms; //!< A shared pointer to the forward transform
		std::vector<std::shared_ptr<solver <datatype> > > solvers; //!< A vector of shared pointers to the matrix solvers
		
		std::vector <std::shared_ptr <plan <datatype> > > pre_transform_plans; //!< A vector of shared pointers of plans to be executed
		std::vector <std::shared_ptr <plan <datatype> > > mid_transform_plans; //!< A vector of shared pointers of plans to be executed
		std::vector <std::shared_ptr <plan <datatype> > > post_transform_plans; //!< A vector of shared pointers of plans to be executed
	};
} /* bases */

#endif /* end of include guard: ELEMENT_HPP_IUTSU4TQ */
