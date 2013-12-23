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
		typedef typename std::map <int, std::shared_ptr<solver <datatype> > >::iterator iterator;
		
		iterator begin () {
			return solvers.begin ();
		}
		
		iterator end () {
			return solvers.end ();
		}
		
		/*!*******************************************************************
		* \param i_name The string representation of the element
		* \param n_boundaries The integer number of boundaries (must be a multiple of 2)
		* \param i_params The parameter object that contains the input parameters of the run
		* \param i_messenger_ptr A pointer to a messenger object
		* \param i_element_flags An integer set of execution element_flags
		*********************************************************************/
		element (int i_name, int i_dimensions, io::parameters <datatype>& i_params, messenger* i_messenger_ptr, int i_element_flags) : 
		params (i_params) {
			element_flags [state] = i_element_flags;
			name = i_name;
			grids.resize (i_dimensions);
			messenger_ptr = i_messenger_ptr;
			timestep = 0.0;
			duration = 0.0;
		}
		
		virtual ~element () {
			// printf ("Destroying bases element\n");
			// printf ("Destroying %p\n", &name);
			// printf ("Destroying %p\n", &messenger_ptr);
			// printf ("Destroying %p\n", &duration);
			// printf ("Destroying %p\n", &timestep);
			// printf ("Destroying %p\n", &scalars);
			// printf ("Destroying %p\n", &grids);
			// printf ("Destroying %p\n", &*failsafe_dump);
			// printf ("Destroying %p\n", &*normal_stream);
			// // printf ("Destroying %p\n", &*transform_stream);
			// printf ("Last\n");
		}
		
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
		 * For simplicity, this could be overloaded in higher dimensions.
		 * 
		 * \return A datatype reference to the given index of the named scalar
		 *********************************************************************/
		virtual datatype& operator() (int name, int index = 0) {
			return (&((*this) [name])) [index];
		}
		
		
		/*!**********************************************************************
		 * \brief Get the pointer to a given named scalar index
		 * 
		 * \param name The integer name from the index enumeration
		 * \param index The integer index of interest
		 * 
		 * For simplicity, this could be overloaded in higher dimensions.
		 * 
		 * \return A pointer to the given index of the named scalar
		 ************************************************************************/
		virtual datatype* ptr (int name, int index = 0) {
			/*
				TODO It would be nice to check if name exists...
			*/
			return &((*this) [name]) + index;
		}
		
		/*!*******************************************************************
		 * \brief Set the collocation grid.
		 * 
		 * \param i_grid A shared_ptr to a grid object
		 * \param index The dimensional index to apply the grid.
		 *********************************************************************/
		inline void set_grid (grid <datatype>* i_grid, int index = 0) {
			grids [index] = std::shared_ptr <grid <datatype>> (i_grid);
		}

		/*!*******************************************************************
		 * \brief Add a matrix solver.
		 * 
		 * \param i_solver A pointer to a solver object
		 *********************************************************************/
		inline void add_solver (int i_name, solver <datatype>* i_solver) {
			solvers [i_name] = std::shared_ptr <solver <datatype>> (i_solver);
			TRACE ("Solver added.");
		}
		
		inline datatype *matrix_ptr (int i_name, int index = 0) {
			/*
				TODO It would be nice to check if name exists...
			*/
			return solvers [i_name]->matrix_ptr (index);
		}

		/*!*******************************************************************
		 * \brief Add a forward horizontal transform operation
		 * 
		 * \param i_plan A pointer to a plan object
		 *********************************************************************/
		inline void add_transform (int i_name, plan <datatype>* i_plan, int i_flags) {
			TRACE ("Adding transform...");
			transforms.push_back (i_name);
			std::shared_ptr <plan <datatype>> plan_ptr = std::shared_ptr <plan <datatype>> (i_plan);
			if (i_flags & forward_horizontal) {
				forward_horizontal_transforms [i_name] = plan_ptr;
			}
			if (i_flags & forward_vertical) {
				forward_vertical_transforms [i_name] = plan_ptr;
			}
			if (i_flags & inverse_horizontal) {
				inverse_horizontal_transforms [i_name] = plan_ptr;
			}
			if (i_flags & inverse_vertical) {
				inverse_vertical_transforms [i_name] = plan_ptr;
			}
		}
		
		/*!*******************************************************************
		 * \brief Initialize the scalar name
		 * 
		 * \param name The integer name index to be initialized
		 * \param initial_conditions The datatype array of initial conditions
		 * \param element_flags A set of binary element_flags
		 *********************************************************************/
		virtual void initialize (int i_name, datatype* initial_conditions = NULL, int i_flags = 0x00) {
			element_flags [i_name] = 0x00;
		}
		
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
			transform (forward_horizontal | forward_vertical);
		}
		
		/*!**********************************************************************
		 * \brief Transform the element forward in the horizontal direction
		 ************************************************************************/
		virtual void transform (int i_flags) {
			TRACE ("Transforming...");
			typedef typename std::vector <int>::iterator iterator; 
			for (iterator iter = transforms.begin (); iter != transforms.end (); iter++) {
				if (i_flags & forward_horizontal) {
					if (!(element_flags [*iter] & transformed_horizontal) && forward_horizontal_transforms [*iter]) {
						forward_horizontal_transforms [*iter]->execute (element_flags [state], element_flags [*iter]);
					}
					element_flags [*iter] |= transformed_horizontal;
				}
				if (i_flags & forward_vertical) {
					if (!(element_flags [*iter] & transformed_vertical) && forward_vertical_transforms [*iter]) {
						forward_vertical_transforms [*iter]->execute (element_flags [state], element_flags [*iter]);
					}
					element_flags [*iter] |= transformed_vertical;
				}
				if (i_flags & inverse_horizontal) {
					if (element_flags [*iter] & transformed_horizontal && inverse_horizontal_transforms [*iter]) {
						inverse_horizontal_transforms [*iter]->execute (element_flags [state], element_flags [*iter]);
					}
					element_flags [*iter] &= ~transformed_horizontal;
				}
				if (i_flags & inverse_vertical) {
					if (element_flags [*iter] & transformed_vertical && inverse_vertical_transforms [*iter]) {
						inverse_vertical_transforms [*iter]->execute (element_flags [state], element_flags [*iter]);
					}
					element_flags [*iter] &= ~transformed_vertical;
				}
			}
		}
		
		/*!**********************************************************************
		 * \brief Output to file
		 ************************************************************************/
		virtual void output () {
			if (normal_stream) {
				normal_stream->to_file ();
			}
			// if (transform_stream) {
			// 	transform_stream->to_file ();
			// }
		}

		/*!**********************************************************************
		 * \brief Factorize all solvers
		 ************************************************************************/
		virtual void factorize () {
			for (iterator iter = begin (); iter != end (); iter++) {
				if (!(element_flags [state] & factorized)) {
					iter->second->factorize ();
				}
				element_flags [iter->first] |= factorized;
			}
		}
		
		/*!**********************************************************************
		 * \brief Execute all solvers
		 ************************************************************************/
		virtual void solve () {
			TRACE ("Beginning solve...");
			datatype t_timestep;
			t_timestep = calculate_timestep ();
			messenger_ptr->min (&t_timestep);
			
			transform (forward_horizontal);

			for (iterator iter = begin (); iter != end (); iter++) {
				iter->second->execute (element_flags [state], element_flags [iter->first]);
			}
			
			duration += timestep;
			INFO ("TOTAL TIME: " << duration);
			if (t_timestep != timestep) {
				element_flags [state] &= ~unchanged_timestep;
				for (iterator iter = begin (); iter != end (); iter++) {
					element_flags [iter->first] &= ~factorized;
				}
				INFO ("Updating timestep: " << t_timestep);
			} else {
				element_flags [state] |= unchanged_timestep;
			}
			timestep = t_timestep;
			transform (forward_vertical);
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
		virtual datatype calculate_timestep () = 0;
		
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
		
		datatype duration; //!< The datatype total simulated time
		datatype timestep; //!< The datatype timestep length

		std::map <int, std::vector <datatype> > scalars; //!< A map of scalar vectors
		flags element_flags;
		std::vector <std::shared_ptr <grid <datatype> > > grids; //!< A vector of shared pointers to the collocation grids
		
		std::shared_ptr <io::output> failsafe_dump; //!< An implementation to dump in case of failure
		std::shared_ptr <io::output> normal_stream; //!< An implementation to output in normal space
		std::shared_ptr <io::output> transform_stream; //!< An implementation to output in transform space
		std::map<int, std::shared_ptr<solver <datatype> > > solvers; //!< A vector of shared pointers to the matrix solvers

	private:
		std::vector <int> transforms;
		
		std::map <int, std::shared_ptr<plan <datatype> > > forward_horizontal_transforms; //!< A vector of shared pointers to the forward horizontal transforms
		std::map <int, std::shared_ptr<plan <datatype> > > inverse_horizontal_transforms; //!< A vector of shared pointers to the inverse horizontal transforms
		std::map <int, std::shared_ptr<plan <datatype> > > forward_vertical_transforms; //!< A vector of shared pointers to the forward vertical transforms
		std::map <int, std::shared_ptr<plan <datatype> > > inverse_vertical_transforms; //!< A vector of shared pointers to the inverse vertical transforms
	};
} /* bases */

#endif /* end of include guard: ELEMENT_HPP_IUTSU4TQ */
