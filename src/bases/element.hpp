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
#include "transform.hpp"
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
		/*!**********************************************************************
		 * \brief Element iterator for iterating through the contained solvers
		 ************************************************************************/
		typedef typename std::vector <int>::iterator iterator;
		
		/*!**********************************************************************
		 * \brief Generate an iterator to iterate through the solvers
		 * 
		 * \return Beginning iterator
		 ************************************************************************/
		iterator begin () {
			return solver_keys.begin ();
		}
		
		/*!**********************************************************************
		 * \brief Generate a finishing iterator for comparison
		 * 
		 * \return Ending iterator
		 ************************************************************************/
		iterator end () {
			return solver_keys.end ();
		}
		
		/*!*******************************************************************
		* \param i_name The string representation of the element
		* \param i_dimensions The integer number of dimensions
		* \param i_params The parameter object that contains the input parameters of the run
		* \param i_messenger_ptr A pointer to a messenger object
		* \param i_element_flags An integer set of execution element_flags
		*********************************************************************/
		element (int i_name, int i_dimensions, io::parameters& i_params, messenger* i_messenger_ptr, int i_element_flags) : 
		dimensions (i_dimensions),
		params (i_params) {
			element_flags [state] = i_element_flags;
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
		 * For simplicity, this could be overloaded in higher dimensions.
		 * 
		 * \return A datatype reference to the given index of the named scalar
		 *********************************************************************/
		virtual datatype& operator() (int name, int index = 0) {
			return (&((*this) [name])) [index];
		}
		
		/*!*******************************************************************
		 * \brief Reset all contained solvers
		 *********************************************************************/
		virtual void explicit_reset () {
			TRACE ("Resetting explicits...");
		}
		
		virtual datatype *_initialize (int i_name, datatype *initial_conditions = NULL, int i_flags = 0x00) = 0;

		/*!*******************************************************************
		 * \brief Initialize the scalar name
		 * 
		 * \param name The integer name index to be initialized
		 * \param initial_conditions The datatype array of initial conditions
		 * \param element_flags A set of binary element_flags
		 *********************************************************************/
		virtual datatype *initialize (int i_name, datatype* initial_conditions = NULL, int i_flags = 0x00) {
			element_flags [i_name] = 0x00;
			for (int i = 0; i < dimensions; ++i) {
				if (!grids [i]) {
					grids [i] = generate_grid (i);
				}
			}
			return _initialize (i_name, initial_conditions, i_flags);
		}
		
		virtual std::shared_ptr <grid <datatype>> generate_grid (int index) = 0;
		
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
		 * \brief Add a matrix solver.
		 * 
		 * \param i_name The integer solver name to add
		 * \param i_solver A pointer to a solver object
		 *********************************************************************/
		inline void add_solver (int i_name, solver <datatype>* i_solver) {
			TRACE ("Adding solver...");
			solvers [i_name] = std::shared_ptr <solver <datatype>> (i_solver);
			solver_keys.push_back (i_name);
			TRACE ("Solver added.");
		}
		
		/*!**********************************************************************
		 * \brief Get the pointer to a matrix from a particular solver
		 * 
		 * \param name The integer solver name from the index enumeration
		 * \param index The integer index of interest in the matrix
		 * 
		 * \return A pointer to the given index of the named solver matrix
		 ************************************************************************/
		inline datatype *matrix_ptr (int i_name, int index = 0) {
			/*
				TODO It would be nice to check if name exists...
			*/
			return solvers [i_name]->matrix_ptr (index);
		}
		
		/*!*******************************************************************
		 * \brief Add a forward horizontal transform operation
		 * 
		 * \param i_name The integer name of the scalar associated with the transform
		 * \param i_plan A pointer to the transform plan
		 * \param i_flags The integer transform flags indicating which transform to associate
		 * 
		 * The flags should come from the transform_flags enumeration.
		 *********************************************************************/
		inline void add_transform (int i_name, master_transform <datatype>* i_transform) {
			TRACE ("Adding transform...");
			bool found = false;
			for (int i = 0; i < (int) transforms.size (); ++i) {
				if (transforms [i] == i_name) {
					found = true;
				}
			}
			if (!found) {
				transforms.push_back (i_name);
			}

			// std::shared_ptr <plan <datatype>> plan_ptr = std::shared_ptr <plan <datatype>> (i_plan);
			// if (i_flags & forward_horizontal) {
			// 	forward_horizontal_transforms [i_name] = plan_ptr;
			// }
			// if (i_flags & forward_vertical) {
			// 	forward_vertical_transforms [i_name] = plan_ptr;
			// }
			// if (i_flags & inverse_horizontal) {
			// 	inverse_horizontal_transforms [i_name] = plan_ptr;
			// }
			// if (i_flags & inverse_vertical) {
			// 	inverse_vertical_transforms [i_name] = plan_ptr;
			// }
			
			master_transforms [i_name] = std::shared_ptr <master_transform <datatype>> (i_transform);
		}
		
		/*!**********************************************************************
		 * \brief Transform the element
		 * 
		 * \param i_flags The integer transform flag indicating the transform direction
		 * 
		 * The flags should come from the transform_flags enumeration.
		 ************************************************************************/
		virtual void transform (int i_flags);
		
		void write_transform_data ();
		void read_transform_data ();
	
		virtual void setup (io::input *input_stream) = 0;
		
		/*!**********************************************************************
		 * \brief Output to file
		 ************************************************************************/
		virtual void output () {
			TRACE ("Writing to file...");
			if (normal_stream) {
				normal_stream->to_file ();
			}
		}

		/*!**********************************************************************
		 * \brief Factorize all solvers
		 ************************************************************************/
		virtual void factorize () {
			TRACE ("Factorizing...");
			for (iterator iter = begin (); iter != end (); iter++) {
				if (!(element_flags [*iter] & factorized)) {
					solvers [*iter]->factorize ();
				}
			}
		}
		
		/*!**********************************************************************
		 * \brief Transform to spectral space and execute all solvers
		 ************************************************************************/
		virtual void solve () {
			TRACE ("Beginning solve...");
			// Execute the solvers
			for (iterator iter = begin (); iter != end (); iter++) {
				solvers [*iter]->execute ();
			}
			// Make certain everything is fully transformed
			write_transform_data ();
			
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
		virtual datatype calculate_min_timestep () = 0;
		
		/*!**********************************************************************
		 * \brief The main function call of the class
		 * 
		 * This method tells the element to begin the main run of the simulation.
		 * It runs through all the specified plans in the appropriate order, and 
		 * updates the values as necessary. Output, if desired, is specified by 
		 * the output streams.
		 ************************************************************************/
		virtual void run (int n_steps);
		
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
		int name, dimensions; //!< An integer representation of the element, to be used in file output
		io::parameters& params; //!< The map that contains the input parameters
		messenger* messenger_ptr; //!< A pointer to the messenger object
		
		datatype duration; //!< The datatype total simulated time
		datatype timestep; //!< The datatype timestep length

		std::map <int, std::vector <datatype> > scalars; //!< A map of scalar vectors
		std::map <int, int> element_flags; //!< A map of integer flags
		std::vector <std::shared_ptr <grid <datatype>>> grids; //!< A vector of shared pointers to the collocation grids
		
		std::shared_ptr <io::output> failsafe_dump; //!< An implementation to dump in case of failure
		std::shared_ptr <io::output> normal_stream; //!< An implementation to output in normal space
		std::shared_ptr <io::output> transform_stream; //!< An implementation to output in transform space
		
		std::vector <int> solver_keys; //!< A vector of integer keys to the solvers map
		std::map<int, std::shared_ptr<solver <datatype>>> solvers; //!< A vector of shared pointers to the matrix solvers

	private:
		std::vector <int> transforms; //!< A vector of integer keys to the transform maps
		std::map <int, std::shared_ptr <master_transform <datatype>>> master_transforms;
	};
} /* bases */

#endif /* end of include guard: ELEMENT_HPP_IUTSU4TQ */
