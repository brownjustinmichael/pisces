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
#include "solver.hpp"
#include "transform.hpp"
#include "../utils/io.hpp"
#include "../utils/messenger.hpp"
#include "collocation.hpp"
#include "../config.hpp"

enum element_flags {
	recv_first = 0x800
};

enum boundary_flags {
	linked_0 = 0x20,
	linked_n = 0x40
};

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
		element (int i_index, int i_name, int n_boundaries, io::parameter_map& i_inputParams, utils::messenger* i_messenger_ptr, int i_flags) : inputParams (i_inputParams) {
			index = i_index;
			name = i_name;
			boundary_bools.resize (n_boundaries);
			boundary_processes.resize (n_boundaries);
			boundary_send_tags.resize (n_boundaries);
			boundary_weights.resize (n_boundaries);
			boundary_recv_tags.resize (n_boundaries);
			boundary_index.resize (n_boundaries);
			excesses.resize (n_boundaries);
			inputParams = i_inputParams;
			messenger_ptr = i_messenger_ptr;
			flags = i_flags;
			timestep = 0.0;
			duration = 0.0;
			logger = config::make_logger (name);
		}
		
		virtual ~element () {
			TRACE (logger, "Calling destructor.");
		}
		
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
		
		int get_index () {
			return index;
		}
		
		bool is_linked (int edge) {
			return boundary_bools [edge];
		}
		
		int get_boundary_index (int edge) {
			return boundary_index [edge];
		}
		
		int get_excess (int edge) {
			return excesses [edge];
		}
		
		int get_expected_excess (int edge) {
			return excesses [edge];
			/*
				TODO Add ability to have different incoming and outgoing excess
			*/
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
				transform ();
			}
		}
		
		virtual void transform () {
			transform_forward->execute ();
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
		inline void set_solver (std::shared_ptr<solver> i_plan) {
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
			TRACE (logger, "Adding plan...");
			plans.push_back (std::move (i_plan));
			TRACE (logger, "Added.");
		}
		
		inline void add_boundary (int edge, int send_tag, int recv_tag, int process, int global_index) {
			TRACE (logger, "Adding boundary, new...");
			boundary_bools [edge] = true;
			/*
				TODO When ready to fix, set boundary_bools [edge] = true;
			*/
			boundary_index [edge] = global_index;
			boundary_weights [edge] = 0.5;
			boundary_send_tags [edge] = send_tag;
			boundary_recv_tags [edge] = recv_tag;
			boundary_processes [edge] = process;
			TRACE (logger, "Added.");
		}
		
		/*!*******************************************************************
		 * \brief Calculate the matrix terms to be used in update
		 * 
		 * In general, this should not be overwritten in subclasses.
		 *********************************************************************/
		virtual void calculate ();
		
		virtual void update_globals (int N, double* global_matrix, double* global_rhs, int* status);
		
		virtual void update_from_globals (double* global_out);
		
		/*!*******************************************************************
		 * \brief Send all relevant boundary data to adjacent elements
		 * 
		 * In general, this should not be overwritten in subclasses.
		 *********************************************************************/
		virtual void send (int n, double* value, int edge, int inc = 1) {
			TRACE (logger, "Sending...");
			if (boundary_bools [edge]) {
				messenger_ptr->send (value, boundary_processes [edge], boundary_send_tags [edge], boundary_weights [edge], n, inc);
			}
		}
		virtual void send (int n, double weight, double* value, int edge, int inc = 1) {
			TRACE (logger, "Sending...");
			if (boundary_bools [edge]) {
				messenger_ptr->send (value, boundary_processes [edge], boundary_send_tags [edge], weight, n, inc);
			}
		}
		virtual void send (int n, int* value, int edge, int inc = 1) {
			TRACE (logger, "Sending...");
			if (boundary_bools [edge]) {
				messenger_ptr->send (value, boundary_processes [edge], boundary_send_tags [edge], boundary_weights [edge], n, inc);
			}
		}
		virtual void send (int n, int weight, int* value, int edge, int inc = 1) {
			TRACE (logger, "Sending...");
			if (boundary_bools [edge]) {
				messenger_ptr->send (value, boundary_processes [edge], boundary_send_tags [edge], weight, n, inc);
			}
		}
		/*!*******************************************************************
		 * \brief Receive all relevant boundary data from adjacent elements
		 * 
		 * In general, this should not be overwritten in subclasses.
		 *********************************************************************/
		virtual void recv (int n, double* value, int edge, int inc = 1) {
			TRACE (logger, "Recving...");
			if (boundary_bools [edge]) {
				messenger_ptr->recv (value, boundary_processes [edge], boundary_recv_tags [edge], boundary_weights [edge], n, inc);
			}
		}
		virtual void recv (int n, double weight, double* value, int edge, int inc = 1) {
			TRACE (logger, "Recving...");
			if (boundary_bools [edge]) {
				messenger_ptr->recv (value, boundary_processes [edge], boundary_recv_tags [edge], weight, n, inc);
			} else {
				*value *= weight;
			}
		}
		virtual void recv (int n, int* value, int edge, int inc = 1) {
			TRACE (logger, "Recving...");
			if (boundary_bools [edge]) {
				messenger_ptr->recv (value, boundary_processes [edge], boundary_recv_tags [edge], boundary_weights [edge], n, inc);
			}
		}
		virtual void recv (int n, int weight, int* value, int edge, int inc = 1) {
			TRACE (logger, "Recving...");
			if (boundary_bools [edge]) {
				messenger_ptr->recv (value, boundary_processes [edge], boundary_recv_tags [edge], weight, n, inc);
			} else {
				*value *= weight;
			}
		}
		
		/*!*******************************************************************
		 * \brief Execute the boundary conditions
		 * 
		 * In general, this should not be overwritten in subclasses.
		 *********************************************************************/
		virtual void execute_boundaries () = 0;
		
		virtual void output ();

		/*!*******************************************************************
		 * \brief Update the element
		 * 
		 * In general, this should not be overwritten in subclasses.
		 *********************************************************************/
		virtual void attempt_update ();
		
		virtual void calculate_bounds ();
		
		virtual void send_positions ();
		
		virtual void recv_positions ();
		
		virtual void send_bounds ();
		
		virtual void recv_bounds ();
		
		virtual void calculate_error ();
		
		virtual void send_error ();
		
		virtual void recv_error ();
		
		virtual void update ();
		
		virtual void update_timestep (double new_timestep);
		
		virtual double calculate_timestep () = 0;
		
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
		int index;
		int name; //!< A string representation of the element, to be used in file output
		io::parameter_map& inputParams; //!< The map that contains the input parameters
		utils::messenger* messenger_ptr;
		
		int flags; //!< An integer set of execution flags
		int logger; //!< An integer representation of the logger, which is interpreted by the config class

		double duration;
		double timestep; //!< The double timestep length

		std::vector <int> names; //!< A vector of integer name indices of the contained scalars
		
		std::shared_ptr<collocation_grid> grid; //!< A shared pointer to the collocation grid
		
		std::shared_ptr<plan> transform_forward; //!< A shared pointer to the forward transform

		std::shared_ptr<io::output> failsafe_dump; //!< An implementation to dump in case of failure
		std::shared_ptr<io::output> normal_stream; //!< An implementation to output in normal space
		std::shared_ptr<io::output> transform_stream; //!< An implementation to output in transform space

		std::vector <bool> boundary_bools;
		std::vector <int> boundary_send_tags;
		std::vector <int> boundary_recv_tags;
		std::vector <int> boundary_processes;
		std::vector <int> boundary_index;
		std::vector <double> boundary_weights;
		std::vector <int> excesses;

	private:
		std::shared_ptr<solver> matrix_solver; //!< A shared pointer to the matrix solver
		
		std::vector <std::shared_ptr <plan>> plans; //!< A vector of shared pointers of plans to be executed
	};
} /* bases */

#endif /* end of include guard: ELEMENT_HPP_IUTSU4TQ */
