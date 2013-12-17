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
			// for (int i = 0; i < (int) forward_horizontal_transforms.size (); ++i) {
			// 	printf ("Destroying %p\n", &*forward_horizontal_transforms [i]);
			// }
			// for (int i = 0; i < (int) inverse_horizontal_transforms.size (); ++i) {
			// 	printf ("Destroying %p\n", &*inverse_horizontal_transforms [i]);
			// }
			// for (int i = 0; i < (int) forward_vertical_transforms.size (); ++i) {
			// 	printf ("Destroying %p\n", &*forward_vertical_transforms [i]);
			// }
			// for (int i = 0; i < (int) inverse_vertical_transforms.size (); ++i) {
			// 	printf ("Destroying %p\n", &*inverse_vertical_transforms [i]);
			// }
			// for (int i = 0; i < (int) solvers.size (); ++i) {
			// 	printf ("Destroying %p\n", &*solvers [i]);
			// }
			// for (int i = 0; i < (int) pre_transform_plans.size (); ++i) {
			// 	printf ("Destroying %p\n", &*pre_transform_plans [i]);
			// }
			// for (int i = 0; i < (int) mid_transform_plans.size (); ++i) {
			// 	printf ("Destroying %p\n", &*mid_transform_plans [i]);
			// }
			// for (int i = 0; i < (int) post_transform_plans.size (); ++i) {
			// 	printf ("Destroying %p\n", &*post_transform_plans [i]);
			// }
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
			return solvers [i_name]->matrix_ptr (index);
		}

		/*!*******************************************************************
		 * \brief Add a forward horizontal transform operation
		 * 
		 * \param i_plan A pointer to a plan object
		 *********************************************************************/
		inline void add_forward_horizontal_transform (plan <datatype>* i_plan) {
			forward_horizontal_transforms.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			// transforms.push_back (i_plan));
		}
		
		/*!*******************************************************************
		 * \brief Add an inverse horizontal transform operation
		 * 
		 * \param i_plan A pointer to a plan object
		 *********************************************************************/
		inline void add_inverse_horizontal_transform (plan <datatype>* i_plan) {
			inverse_horizontal_transforms.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			// transforms.push_back (i_plan));
		}
		
		/*!*******************************************************************
		 * \brief Add a forward vertical transform operation
		 * 
		 * \param i_plan A pointer to a plan object
		 *********************************************************************/
		inline void add_forward_vertical_transform (plan <datatype>* i_plan) {
			forward_vertical_transforms.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			// transforms.push_back (i_plan));
		}
		
		/*!*******************************************************************
		 * \brief Add an inverse vertical transform operation
		 * 
		 * \param i_plan A pointer to a plan object
		 *********************************************************************/
		inline void add_inverse_vertical_transform (plan <datatype>* i_plan) {
			inverse_vertical_transforms.push_back (std::shared_ptr <plan <datatype>> (i_plan));
			// transforms.push_back (i_plan));
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
		
		/*!*******************************************************************
		 * \brief Initialize the scalar name
		 * 
		 * \param name The integer name index to be initialized
		 * \param initial_conditions The datatype array of initial conditions
		 * \param flags A set of binary flags
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
		 * \brief Transform the element forward in the horizontal direction
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
		
		/*!**********************************************************************
		 * \brief Transform the element backward in the horizontal direction
		 ************************************************************************/
		virtual void transform_horizontal_inverse () {
			TRACE ("Inverting...");
			if (flags & transformed_horizontal) {
				for (int i = 0; i < (int) inverse_horizontal_transforms.size (); ++i) {
					DEBUG ("PrePoint " << i);
					DEBUG ("Point " << &*(inverse_horizontal_transforms [i]));
					inverse_horizontal_transforms [i]->execute (flags);
				}
				flags &= ~transformed_horizontal;
			}
			TRACE ("Done.");
		}
		
		/*!**********************************************************************
		 * \brief Transform the element forward in the vertical direction
		 ************************************************************************/
		virtual void transform_vertical_forward () {
			TRACE ("Transforming...");
			if (!(flags & transformed_vertical)) {
				for (int i = 0; i < (int) forward_vertical_transforms.size (); ++i) {
					forward_vertical_transforms [i]->execute (flags);
				}
				flags |= transformed_vertical;
			}
		}
		
		/*!**********************************************************************
		 * \brief Transform the element backward in the vertical direction
		 ************************************************************************/
		virtual void transform_vertical_inverse () {
			TRACE ("Inverting...");
			if (flags & transformed_vertical) {
				for (int i = 0; i < (int) inverse_vertical_transforms.size (); ++i) {
					inverse_vertical_transforms [i]->execute (flags);
				}
				flags &= ~transformed_vertical;
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
			if (!(flags & factorized)) {
				typedef typename std::map<int, std::shared_ptr<solver <datatype> > >::iterator iterator; 
				for (iterator iter = solvers.begin (); iter != solvers.end (); iter++) {
					iter->second->factorize ();
				}
				flags |= factorized;
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
			
			typedef typename std::map<int, std::shared_ptr<solver <datatype> > >::iterator iterator; 
			for (iterator iter = solvers.begin (); iter != solvers.end (); iter++) {
				iter->second->execute (flags);
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

		std::map <int, std::vector <datatype> > scalars; //!< A map of scalar vectors
		std::vector <std::shared_ptr <grid <datatype> > > grids; //!< A vector of shared pointers to the collocation grids
		
		std::shared_ptr <io::output> failsafe_dump; //!< An implementation to dump in case of failure
		std::shared_ptr <io::output> normal_stream; //!< An implementation to output in normal space
		std::shared_ptr <io::output> transform_stream; //!< An implementation to output in transform space

	private:
		std::vector<std::shared_ptr<plan <datatype> > > forward_horizontal_transforms; //!< A vector of shared pointers to the forward horizontal transforms
		std::vector<std::shared_ptr<plan <datatype> > > inverse_horizontal_transforms; //!< A vector of shared pointers to the inverse horizontal transforms
		std::vector<std::shared_ptr<plan <datatype> > > forward_vertical_transforms; //!< A vector of shared pointers to the forward vertical transforms
		std::vector<std::shared_ptr<plan <datatype> > > inverse_vertical_transforms; //!< A vector of shared pointers to the inverse vertical transforms
		std::map<int, std::shared_ptr<solver <datatype> > > solvers; //!< A vector of shared pointers to the matrix solvers
		
		/*
			TODO Make solvers a map so that matrices are easy access
		*/
		
		std::vector <std::shared_ptr <plan <datatype> > > pre_transform_plans; //!< A vector of shared pointers of plans to be executed before the transforms
		std::vector <std::shared_ptr <plan <datatype> > > mid_transform_plans; //!< A vector of shared pointers of plans to be executed after the vertical transform
		std::vector <std::shared_ptr <plan <datatype> > > post_transform_plans; //!< A vector of shared pointers of plans to be executed after both transforms
	};
} /* bases */

#endif /* end of include guard: ELEMENT_HPP_IUTSU4TQ */
