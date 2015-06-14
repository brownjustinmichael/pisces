/*!***********************************************************************
 * \file element.hpp
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

#include "mpi/messenger.hpp"

#include <string>
#include <cassert>
#include <memory>
#include <pthread.h>

#include <gsl/gsl_siman.h>
#include <omp.h>

#include "logger/logger.hpp"
#include "versions/version.hpp"
#include "io/input.hpp"
#include "io/output.hpp"
#include "io/parameters.hpp"
#include "io/formats/virtual.hpp"
#include "plans/grids/grid.hpp"
#include "plans/plan.hpp"
#include "plans-solvers/equation.hpp"
#include "plans-transforms/transformer.hpp"

#include "rezone.hpp"
#include "data/data.hpp"

/*!**********************************************************************
 * \namespace pisces
 * 
 * \brief A namespace containing the element class, the main class of the code
 ************************************************************************/
namespace pisces
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
	protected:
		std::vector <grids::axis> axes; //!< A vector of axis objects, containing the basic grid information
		int name; //!< An integer representation of the element, to be used in file output
		int dimensions; //!< The integer number of dimensions in the element
		io::parameters& params; //!< The map that contains the input parameters
		
		datatype &timestep; //!< The datatype timestep length

		data::data <datatype> &data; //!< An object that contains all the data in the simulation
		std::map <std::string, int> &element_flags; //!< A map of integer flags
		
		std::map <std::string, std::shared_ptr <plans::solvers::equation <datatype>>> equations; //!< A vector of shared pointers to the matrix equations
		std::map <std::string, std::shared_ptr <plans::transforms::transformer <datatype>>> transformers; //!< A map containing the transformer objects for each variable
		std::vector <std::string> transforms; //!< A vector containing the variables with transforms
		int transform_threads; //!< The number of transform threads available
		
	private:
		datatype rezone_mult; //!< To merit rezoning, the new timestep must be at least this factor larger than the current one
		std::vector <std::string> equation_keys; //!< A vector of integer keys to the equations map
		
	public:
		datatype &duration; //!< The datatype total simulated time
		std::vector <std::shared_ptr <grids::grid <datatype>>> grids; //!< A vector of shared pointers to the collocation grids
		mpi::messenger* messenger_ptr; //!< A pointer to the messenger object
		formats::virtual_file *rezone_virtual_file; //!< A shared_ptr to a virtual file object, for rezoning
		/*!**********************************************************************
		 * \brief Element iterator for iterating through the contained equations
		 ************************************************************************/
		typedef typename std::vector <std::string>::iterator iterator;
		
		/*!**********************************************************************
		 * \brief Generate an iterator to iterate through the equations
		 * 
		 * \return Beginning iterator
		 ************************************************************************/
		iterator begin () {
			return equation_keys.begin ();
		}
		
		/*!**********************************************************************
		 * \brief Generate a finishing iterator for comparison
		 * 
		 * \return Ending iterator
		 ************************************************************************/
		iterator end () {
			return equation_keys.end ();
		}
		
		/*!*******************************************************************
		* \param i_name The string representation of the element
		* \param i_dimensions The integer number of dimensions
		* \param i_params The parameter object that contains the input parameters of the run
		* \param i_data An object that contains all the data in the simulation
		* \param i_messenger_ptr A pointer to a messenger object for inter-element communication
		* \param i_element_flags An integer set of global flags for the element
		*********************************************************************/
		element (int i_name, int i_dimensions, io::parameters& i_params, data::data <datatype> &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) : 
		dimensions (i_dimensions),
		params (i_params),
		timestep (i_data.timestep),
		data (i_data),
		element_flags (data.flags),
		duration (i_data.duration) {
			element_flags ["element"] = i_element_flags;
			name = i_name;
			grids.resize (i_dimensions);
			axes.resize (i_dimensions);
			messenger_ptr = i_messenger_ptr;
			timestep = 0.0;
			transform_threads = params.get <int> ("parallel.transform.threads");
			rezone_mult = params.get <datatype> ("grid.rezone.mult");
		}
		
		virtual ~element () {}
		
		/*!**********************************************************************
		 * \brief Gets the version of the class
		 * 
		 * \return The version of the class
		 ************************************************************************/
		static versions::version& version () {
			static versions::version version ("1.1.0.0");
			return version;
		}
		
		/*!*******************************************************************
		 * \brief Get the datatype reference to the named scalar
		 * 
		 * \param name The integer name from the index enumeration
		 * 
		 * \return A datatype reference to the first element of the named scalar
		 *********************************************************************/
		inline datatype *operator[] (std::string name) {
			// if (scalars.find (name) == scalars.end ()) {
			// 	FATAL ("Index " << name << " not found in element.");
			// 	throw 0;
			// }
			// return scalars [name] [0];
			DEBUG (data [name]);
			return data [name];
		}
		
		/*!*******************************************************************
		 * \brief Get the datatype reference to the given index of the named scalar
		 * 
		 * \param name The integer name from the index enumeration
		 * \param index The integer index of interest	
		 * 
		 * For simplicity, this could be overwritten in higher dimensions.
		 * 
		 * \return A datatype reference to the given index of the named scalar
		 *********************************************************************/
		virtual datatype& operator() (std::string name, int index = 0) {
			return ((*this) [name]) [index];
		}
		
		/*!**********************************************************************
		 * \brief Get the mode integer for the given element
		 * 
		 * This must be overwritten in elements. It is to make certain that loaded inputs have the same spatial geometry as the element that they are being loaded into.
		 ************************************************************************/
		virtual const int &get_mode () = 0;
		
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
		virtual datatype* ptr (std::string name, int index = 0) {
			/*
				TODO It would be nice to check if name exists in debug mode...
			*/
			// return &((*this) [name]) + index;
			return data (name, index);
		}
		
		/*!*******************************************************************
		 * \brief Add a matrix equation.
		 * 
		 * \param i_name The integer solver name to add
		 * \param i_solver_ptr A pointer to a solver object
		 *********************************************************************/
		inline void add_equation (std::string i_name, std::shared_ptr <plans::solvers::equation <datatype> > i_solver_ptr) {
			TRACE ("Adding solver...");
			equations [i_name] = i_solver_ptr;
			equation_keys.push_back (i_name);
			TRACE ("Solver added.");
		}
		
		/*!**********************************************************************
		 * \brief Transform the element according to the flags
		 * 
		 * \param i_flags The flags to pass to the transformer transform member (e.g. forward_horizontal, read_before)
		 ************************************************************************/
		void transform (int i_flags) {
			TRACE ("Transforming...");
			int threads = transform_threads;
			#pragma omp parallel for num_threads (threads)
			for (int i = 0; i < (int) transforms.size (); ++i) {
				transformers [transforms [i]]->transform (i_flags);
			}
		}
	
		/*!**********************************************************************
		 * \brief Get the pointer to a matrix from a particular solver
		 * 
		 * \param i_name The integer solver name from the index enumeration
		 * \param index The integer index of interest in the matrix
		 * 
		 * \return A pointer to the given index of the named solver matrix
		 ************************************************************************/
		inline datatype *matrix_ptr (std::string i_name, int index = 0) {
			/*
				TODO It would be nice to check if name exists in debug mode...
			*/
			return equations [i_name]->matrix_ptr (index);
		}
		
		/*!*******************************************************************
		 * \brief Initialize the scalar name and generate its transforms
		 * 
		 * \param i_name The integer index to be initialized from the index enumeration
		 *********************************************************************/
		virtual datatype *initialize (std::string i_name) {
			return this->ptr (i_name);
		}
		
		/*!**********************************************************************
		 * \brief Generate a shared_ptr to a grid for a specified dimension (use last dimension if unspecified)
		 * 
		 * \param axis_ptr A pointer to an axis object, which contains the extent of the grid and number of gridpoints
		 * \param index The index of the grid to add (defaults to the next available)
		 ************************************************************************/
		virtual std::shared_ptr <grids::grid <datatype>> generate_grid (grids::axis *axis_ptr, int index = -1) = 0;
		
		/*!**********************************************************************
		 * \brief Factorize all equations
		 ************************************************************************/
		virtual void factorize () {
			TRACE ("Factorizing...");
			for (iterator iter = begin (); iter != end (); iter++) {
				if (!(element_flags [*iter] & plans::solvers::factorized)) {
					// DEBUG ("Factorizing " << *iter);
					equations [*iter]->factorize ();
				}
			}
		}
		
		/*!**********************************************************************
		 * \brief Transform to spectral space and execute all equations
		 ************************************************************************/
		virtual void solve () {
			TRACE ("Beginning solve...");
			// Execute the equations
			// std::map <int, omp_lock_t> locks;
			for (iterator iter = begin (); iter != end (); iter++) {
				element_flags [*iter] &= ~solved;
			}
			bool completely_solved = false, skip;
			std::vector <std::string> can_be_solved;
			
			while (!completely_solved) {
				can_be_solved.clear ();
				for (iterator iter = begin (); iter != end (); iter++) {
					if (element_flags [*iter] & solved) {
						continue;
					}
					
					skip = false;
					for (int j = 0; j < equations [*iter]->n_dependencies (); ++j) {
						if (!(*(equations [*iter]->get_dependency (j)->component_flags) & solved)) {
							skip = true;
							break;
						}
					}
					
					if (!skip) {
						can_be_solved.push_back (*iter);
					}
				}
				
				int threads = std::min (params.get <int> ("parallel.solver.threads"), (int) can_be_solved.size ());
				// #pragma omp parallel for num_threads (threads)
				for (int i = 0; i < threads; ++i) {
					std::string name = can_be_solved [i];
					equations [name]->solve ();
					// #pragma omp atomic
					element_flags [name] |= solved;
				}
				
				completely_solved = true;
				for (iterator iter = begin (); iter != end (); iter++) {
					completely_solved = completely_solved && (element_flags [*iter] & solved);
				}
			}
			
			for (iterator iter = begin (); iter != end (); iter++) {
				equations [*iter]->reset ();
			}
			
			// Make certain everything is fully transformed
			
			
			transform (plans::transforms::forward_vertical | plans::transforms::no_read);
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
		virtual datatype calculate_min_timestep (formats::virtual_file *virtual_file = NULL, bool limiters = true) = 0;
		
		/*!**********************************************************************
		 * \brief Caculate the zoning that minimizes the timestep according to size constraints
		 * 
		 * \param positions A datatype pointer to the zoning positions
		 * \param min_size The datatype minimum distance between zoning positions
		 * \param max_size The datatype maximum distance between zoning positions
		 * \param n_tries the number of attempts before stepping
		 * \param iters_fixed_t the integer number of iterations at each temperature
		 * \param step_size the real maximum step size
		 * \param k a real boltzmann constant
		 * \param t_initial the real initial temperature
		 * \param mu_t the real damping factor for temperature
		 * \param t_min the real damping factor parameter for temperature
		 * 
		 * Using a simulated annealing technique, rezone all the elements such that the timestep is minimized across them. This is an expensive operation and should be used sparingly.
		 * 
		 * \return A shared_ptr to a virtual_file object containing the chosen rezoning
		 ************************************************************************/
		virtual formats::virtual_file *rezone_minimize_ts (datatype *positions, datatype min_size, datatype max_size, int n_tries = 20, int iters_fixed_t = 1000, datatype step_size = 1.0, datatype k = 1.0, datatype t_initial = 0.008, datatype mu_t = 1.003, datatype t_min = 2.0e-6, double (*function) (element <datatype> *, formats::virtual_file *) = element <datatype>::rezone_calculate_ts) {
			TRACE ("Rezoning...");
			transform (plans::transforms::inverse_horizontal | plans::transforms::inverse_vertical);

			rezone_virtual_file = make_virtual_file (profile_only | timestep_only);
			
			rezone_data <datatype> data;
			data.element_ptr = this;
			data.id = messenger_ptr->get_id ();
			data.np = messenger_ptr->get_np ();
			data.merit_func = function;
			data.min_size = min_size;
			data.max_size = max_size;
			
			get_zoning_positions (data.positions);
			
			datatype original = rezone_data <datatype>::func (&data);
			
			gsl_siman_params_t params = {n_tries, iters_fixed_t, step_size, k, t_initial, mu_t, t_min};

			const gsl_rng_type * T;
			gsl_rng * r;
			
			gsl_rng_env_setup();
			
			T = gsl_rng_default;
			r = gsl_rng_alloc(T);

			gsl_siman_solve(r, &data, rezone_data <datatype>::func, rezone_data <datatype>::rezone_generate_step, rezone_data <datatype>::rezone_step_size, NULL, NULL, NULL, NULL, sizeof (rezone_data <datatype>), params);

			gsl_rng_free (r);
			
			datatype value = rezone_data <datatype>::func (&data);
			
			DEBUG ("Old: " << original << " New: " << value);
			
			if (value < original * (original < 0. ? 1. / rezone_mult : rezone_mult)) {
				for (int i = 0; i < data.np + 1; ++i) {
					positions [i] = data.positions [i];
				}
				return make_rezoned_virtual_file (data.positions, make_virtual_file ());
			}
			
			return NULL;
		}
		
		/*!**********************************************************************
		 * \brief The main function call of the class
		 * 
		 * This method tells the element to begin the main run of the simulation. It runs through all the specified plans in the appropriate order, and updates the values as necessary. Output, if desired, is specified by the output streams.
		 ************************************************************************/
		virtual void run (int &n_steps, int max_steps, int check_every = -1, datatype stop = 100.0);
		
		/*!**********************************************************************
		 * \brief Protected method to make a virtual file of the current state
		 * 
		 * \param flags The binary flags to pass to the method
		 * 
		 * This method needs to be overwritten in a subclass.
		 * 
		 * \return A shared_ptr to the virtual_file of the current state
		 ************************************************************************/
		virtual formats::virtual_file *make_virtual_file (int flags = 0x00) = 0;
	
		/*!**********************************************************************
		 * \brief Protected method to rezone the current state into a virtual file
		 * 
		 * \param positions A datatype pointer to the zoning position array
		 * \param virtual_file_ptr A pointer to the virtual_file to be rezoned (for speed)
		 * \param flags The binary flags for method execution
		 * 
		 * The method needs to be overwritten in a subclass
		 * 
		 * \return A shared_ptr to the rezoned virtual file
		 ************************************************************************/
		virtual formats::virtual_file *make_rezoned_virtual_file (datatype *positions, formats::virtual_file *virtual_file_ptr, int flags = 0x00) = 0;
		
	protected:
		/*!**********************************************************************
		 * \brief Get the current zoning position array
		 * 
		 * \param positions A datatype pointer to the input and output position array
		 * 
		 * This method needs to be overwritten in a subclass
		 ************************************************************************/
		virtual void get_zoning_positions (datatype *positions) = 0;
		
	private:
		/*!**********************************************************************
		 * \brief Given a set of zoning data, calculate the timestep
		 * 
		 * \param i_rezone_data A pointer to an array of rezone_union objects
		 * 
		 * In order, the rezone_union objects must contain 1) a pointer to the element, 2) the total number of elements in the system, 3) min_size argument from rezone_minimize_ts, 4) max_size argument from rezone_minimize_ts, 5-n) positions of the zonal boundaries. This method is to be called by the GSL simulated annealing routine.
		 ************************************************************************/
		static double rezone_calculate_ts (element <datatype> *element_ptr, formats::virtual_file *virt) {
			double timestep = element_ptr->calculate_min_timestep (virt, false);
			return -timestep;
		}
		
		/*
			TODO Allow for upside-down elements
		*/
	};
} /* pisces */

#endif /* end of include guard: ELEMENT_HPP_IUTSU4TQ */
