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

#include "../utils/messenger.hpp"
#include <string>
#include <cassert>
#include <memory>
#include "plan.hpp"
#include "solver.hpp"
#include "../utils/io.hpp"
#include "grid.hpp"
#include "transform.hpp"
#include "../config.hpp"
#include <gsl/gsl_siman.h>

namespace bases
{
	template <class datatype>
	class element;
	
	/*!**********************************************************************
	 * \brief This is a simple union for the zoning minimization
	 ************************************************************************/
	template <class datatype>
	union rezone_union {
		bases::element <datatype> *element_ptr; //! A pointer to an element
		int np; //! The integer number of elements
		datatype position; //! A real zoning position
	};
	
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
		std::vector <bases::axis> axes; //!< A vector of axis objects, containing the basic grid information
		utils::messenger* messenger_ptr; //!< A pointer to the messenger object
		int name; //!< An integer representation of the element, to be used in file output
		int dimensions; //!< The integer number of dimensions in the element
		io::parameters& params; //!< The map that contains the input parameters
		
		datatype duration; //!< The datatype total simulated time
		datatype timestep; //!< The datatype timestep length

		std::vector <std::shared_ptr <grid <datatype>>> grids; //!< A vector of shared pointers to the collocation grids
		std::map <int, std::vector <datatype> > scalars; //!< A map of scalar vectors
		std::map <int, std::string> scalar_names; //!< A map of string representations of the scalars
		std::map <int, int> element_flags; //!< A map of integer flags
		
		std::map<int, std::shared_ptr <master_solver <datatype>>> solvers; //!< A vector of shared pointers to the matrix solvers
		
	private:
		std::shared_ptr <io::output> normal_stream; //!< An implementation to output in normal space
		std::shared_ptr <io::output> transform_stream; //!< An implementation to output in transform space
		std::vector <std::shared_ptr <io::output>> normal_profiles; //!< An implementation to output in transform space

		std::vector <int> solver_keys; //!< A vector of integer keys to the solvers map
		io::virtual_dump *rezone_dump; //!< A shared_ptr to a virtual dump object, for rezoning
		
	public:
		std::vector <int> transforms; //!< A vector of integer keys to the transform maps
		std::map <int, std::shared_ptr <master_transform <datatype>>> master_transforms; //!< A map of shared_ptrs to the transform objects
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
		* \param i_messenger_ptr A pointer to a messenger object for inter-element communication
		* \param i_element_flags An integer set of global flags for the element
		*********************************************************************/
		element (int i_name, int i_dimensions, io::parameters& i_params, utils::messenger* i_messenger_ptr, int i_element_flags) : 
		dimensions (i_dimensions),
		params (i_params) {
			element_flags [state] = i_element_flags;
			name = i_name;
			grids.resize (i_dimensions);
			axes.resize (i_dimensions);
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
		 * For simplicity, this could be overwritten in higher dimensions.
		 * 
		 * \return A datatype reference to the given index of the named scalar
		 *********************************************************************/
		virtual datatype& operator() (int name, int index = 0) {
			return (&((*this) [name])) [index];
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
		virtual datatype* ptr (int name, int index = 0) {
			/*
				TODO It would be nice to check if name exists in debug mode...
			*/
			return &((*this) [name]) + index;
		}
		
		/*!*******************************************************************
		 * \brief Add a matrix solver.
		 * 
		 * \param i_name The integer solver name to add
		 * \param i_solver_ptr A pointer to a solver object
		 *********************************************************************/
		inline void add_solver (int i_name, std::shared_ptr <master_solver <datatype> > i_solver_ptr) {
			TRACE ("Adding solver...");
			solvers [i_name] = i_solver_ptr;
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
				TODO It would be nice to check if name exists in debug mode...
			*/
			return solvers [i_name]->matrix_ptr (index);
		}
		
		/*!*******************************************************************
		 * \brief Initialize the scalar name and generate its transforms
		 * 
		 * \param name The integer index to be initialized from the index enumeration
		 * \param i_str The string representation of the scalar field, for output
		 * \param initial_conditions The datatype array of initial conditions
		 * \param element_flags A set of binary element_flags for instantiation
		 *********************************************************************/
		virtual datatype *initialize (int i_name, std::string i_str, datatype* initial_conditions = NULL, int i_flags = 0x00) {
			/*
				TODO List possible flags in all documentation
			*/
			element_flags [i_name] = 0x00;
			for (int i = 0; i < dimensions; ++i) {
				if (!grids [i]) {
					grids [i] = generate_grid (&axes [i], i);
				}
			}
			scalar_names [i_name] = i_str;
			return _initialize (i_name, initial_conditions, i_flags);
		}
		
		/*!**********************************************************************
		 * \brief Generate a shared_ptr to a grid for a specified dimension (use last dimension if unspecified)
		 * 
		 * \param axis_ptr A pointer to an axis object, which contains the extent of the grid and number of gridpoints
		 * \param 
		 ************************************************************************/
		virtual std::shared_ptr <grid <datatype>> generate_grid (bases::axis *axis_ptr, int index = -1) = 0;
		
		/*!**********************************************************************
		 * \brief Given an input stream, load all the relevant data into the element
		 * 
		 * \param input_stream_ptr A pointer to an input object
		 ************************************************************************/
		void setup (io::input *input_ptr) {
			// Iterate through the scalar fields and append them to the variables for which the input will search
			typedef typename std::map <int, std::vector <datatype> >::iterator iterator;
			for (iterator iter = scalars.begin (); iter != scalars.end (); ++iter) {
				input_ptr->template append <datatype> (scalar_names [iter->first], ptr (iter->first));
			}
			// Also look for the scalars for total simulated time and geometry mode
			input_ptr->template append_scalar <datatype> ("t", &duration);
			int mode;
			input_ptr->template append_scalar <int> ("mode", &mode);
			
			try {
				// Read from the input into the element variables
				input_ptr->from_file ();
			} catch (exceptions::io::bad_variables &except) {
				// Send a warning if there are missing variables in the input
				WARN (except.what ());
			}
			
			// Make sure that the mode matched the mode of the element
			if (mode != get_mode ()) {
				FATAL ("Loading simulation in different mode: " << mode << " instead of " << get_mode ());
				throw 0;
			}
		}
		
		/*!**********************************************************************
		 * \brief Given an output_stream, append all relevant outputs
		 * 
		 * \param output_ptr A shared_ptr to the output object
		 * \param flags The integer flags to specify the type of output (transform stream, normal_stream)
		 * 
		 * Given an output object, prepare it to output all the scalar fields tracked directly by the element. Make it one of the output streams called during output. If flags is normal_stream, output the variables when in Cartesian space, else if flags is transform_stream, output with horizontal grid in Fourier space, but vertical grid in Cartesian space.
		 ************************************************************************/
		void setup_output (std::shared_ptr <io::output> output_ptr, int flags = 0x00) {
			// Iterate through the scalar fields and append them to the variables which the output will write to file
			typedef typename std::map <int, std::vector <datatype> >::iterator iterator;
			for (iterator iter = scalars.begin (); iter != scalars.end (); ++iter) {
				output_ptr->template append <datatype> (scalar_names [iter->first], ptr (iter->first));
			}
			// Also prepare to output the total simulated time and geometry mode
			output_ptr->template append_scalar <datatype> ("t", &duration);
			output_ptr->template append_scalar <const int> ("mode", &(get_mode ()));

			// Check the desired output time and save the output object in the appropriate variable
			if (flags & transform_output) {
				transform_stream = output_ptr;
			} else if (flags & normal_output) {
				normal_stream = output_ptr;
			}
		}
		
		/*!**********************************************************************
		 * \brief Given an output stream, append all relevant outputs for profiling
		 * 
		 * \param output_ptr A shared_ptr to the output object
		 * \param flags The integer flags to specify the type of output profile
		 * 
		 * Given an output object, prepare it to output all the profiles of the scalar fields tracked directly by the element. Make it one of the output streams called during output. This method must be overwritten in subclasses and should be concluded by calling this method.
		 ************************************************************************/
		virtual void setup_profile (std::shared_ptr <io::output> output_ptr, int flags = 0x00) {
			normal_profiles.push_back (output_ptr);
		}
		
		/*
			TODO These three methods are related in purpose but varied in design; they should be unified
		*/
		
		/*!**********************************************************************
		 * \brief Output to file
		 ************************************************************************/
		virtual void output () {
			TRACE ("Writing to file...");
			if (normal_stream) {
				normal_stream->to_file ();
			}
		}
		
		/*!*******************************************************************
		 * \brief Add a transform operation to a specific scalar fields
		 * 
		 * \param i_name The integer name of the scalar associated with the transform
		 * \param i_transform_ptr A shared_ptr to the master_transform plan
		 * 
		 * The master_transform class contains both horizontal and vertical transforms in both directions. It can also send all the transformations to the GPU if desired.
		 *********************************************************************/
		inline void add_transform (int i_name, std::shared_ptr <master_transform <datatype>> i_transform_ptr) {
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
			
			master_transforms [i_name] = i_transform_ptr;
		}
		
		/*!**********************************************************************
		 * \brief Transform the element
		 * 
		 * \param i_flags The integer transform flag indicating the transform direction
		 * 
		 * The flags should come from the transform_flags enumeration.
		 ************************************************************************/
		virtual void transform (int i_flags);
		
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
				TRACE ("Solving " << *iter);
				try {
					solvers [*iter]->solve ();
				} catch (std::exception &except) {
					FATAL ("Failure in solver " << *iter);
					throw except;
				}
			}
			// Make certain everything is fully transformed
			
			transform (forward_vertical | no_read);
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
		virtual datatype calculate_min_timestep (io::virtual_dump *dump = NULL) = 0;
		
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
		 * \return A shared_ptr to a virtual_dump object containing the chosen rezoning
		 ************************************************************************/
		virtual io::virtual_dump *rezone_minimize_ts (datatype * positions, datatype min_size, datatype max_size, int n_tries = 20, int iters_fixed_t = 1000, datatype step_size = 1.0, datatype k = 1.0, datatype t_initial = 0.008, datatype mu_t = 1.003, datatype t_min = 2.0e-6) {
			TRACE ("Rezoning...");
			transform (inverse_horizontal | inverse_vertical);

			rezone_dump = make_dump (profile_only | timestep_only);
			
			rezone_union <datatype> rezone_data [(messenger_ptr->get_np () + 5)];
			rezone_data [0].element_ptr = this;
			rezone_data [1].np = messenger_ptr->get_np ();
			rezone_data [2].position = min_size;
			rezone_data [3].position = max_size;
			
			get_zoning_positions (positions);
			
			for (int i = 0; i < messenger_ptr->get_np () + 1; ++i) {
				rezone_data [i + 4].position = positions [i];
			}
			
			gsl_siman_params_t params = {n_tries, iters_fixed_t, step_size, k, t_initial, mu_t, t_min};

		    const gsl_rng_type * T;
		    gsl_rng * r;

		    gsl_rng_env_setup();

		    T = gsl_rng_default;
		    r = gsl_rng_alloc(T);
			
		    gsl_siman_solve(r, rezone_data, element <datatype>::rezone_calculate_ts, element <datatype>::rezone_generate_step, element <datatype>::rezone_step_size, NULL, NULL, NULL, NULL, sizeof(rezone_union <datatype>) * (messenger_ptr->get_np () + 5), params);

		    gsl_rng_free (r);
			
			if (rezone_calculate_ts (rezone_data) < -timestep) {
				for (int i = 0; i < messenger_ptr->get_np () + 1; ++i) {
					positions [i] = rezone_data [i + 4].position;
				}
			}
			
			return make_rezoned_dump (positions, make_dump ());
		}
		
		/*!**********************************************************************
		 * \brief The main function call of the class
		 * 
		 * This method tells the element to begin the main run of the simulation. It runs through all the specified plans in the appropriate order, and updates the values as necessary. Output, if desired, is specified by the output streams.
		 ************************************************************************/
		virtual void run (int &n_steps, int max_steps, int check_every = -1);
		
	protected:
		/*!**********************************************************************
		 * \brief Protected method to initialize a scalar field and generate its transforms
		 * 
		 * \param i_name The integer reference to the scalar field
		 * \param initial_conditions Any initial conditions to instantiate with, NULL to use default
		 * \param i_flags Integer flags describing the instantiation (e.g. whether the field doesn't require transforms)
		 * 
		 * This method is to be overwritten in subclasses to initialize a scalar field in accordance with the geometry of the element.
		 * 
		 * \return Returns a datatype pointer to the initialized scalar
		 ************************************************************************/
		virtual datatype *_initialize (int i_name, datatype *initial_conditions = NULL, int i_flags = 0x00) = 0;
		
		/*!**********************************************************************
		 * \brief Protected method to make a virtual dump of the current state
		 * 
		 * \param flags The binary flags to pass to the method
		 * 
		 * This method needs to be overwritten in a subclass.
		 * 
		 * \return A shared_ptr to the dump of the current state
		 ************************************************************************/
		virtual io::virtual_dump *make_dump (int flags = 0x00) = 0;
		
		/*!**********************************************************************
		 * \brief Protected method to rezone the current state into a virtual dump
		 * 
		 * \param positions A datatype pointer to the zoning position array
		 * \param dump_ptr A pointer to the dump to be rezoned (for speed)
		 * \param flags The binary flags for method execution
		 * 
		 * The method needs to be overwritten in a subclass
		 * 
		 * \return A shared_ptr to the rezoned virtual dump
		 ************************************************************************/
		virtual io::virtual_dump *make_rezoned_dump (datatype *positions, io::virtual_dump *dump_ptr, int flags = 0x00) = 0;
		
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
		static double rezone_calculate_ts (void *i_rezone_data) {
			bases::element <datatype> *element_ptr = ((rezone_union <datatype> *) i_rezone_data)->element_ptr;
			datatype positions [((rezone_union <datatype> *) i_rezone_data) [1].np + 1];
			utils::messenger *messenger_ptr = element_ptr->messenger_ptr;
			for (int i = 0; i < ((rezone_union <datatype> *) i_rezone_data) [1].np + 1; ++i) {
				positions [i] = ((rezone_union <datatype> *) i_rezone_data) [i + 4].position;
			}
			double timestep = element_ptr->calculate_min_timestep (element_ptr->make_rezoned_dump (positions, &*(element_ptr->rezone_dump), profile_only));
			messenger_ptr->min (&timestep);
			return -timestep;
		}
	
		/*!**********************************************************************
		 * \brief Given a set of zoning data, print the state, for debugging
		 * 
		 * \param i_rezone_data A pointer to an array of rezone_union objects
		 * 
		 * In order, the rezone_union objects must contain 1) a pointer to the element, 2) the total number of elements in the system, 3) min_size argument from rezone_minimize_ts, 4) max_size argument from rezone_minimize_ts, 5-n) positions of the zonal boundaries. This method is to be called by the GSL simulated annealing routine.
		 ************************************************************************/
		static void print_rezone (void *i_rezone_data) {
			for (int i = 0; i < ((rezone_union <datatype> *) i_rezone_data) [1].np + 1; ++i) {
				printf (" %f ", ((rezone_union <datatype> *) i_rezone_data) [i + 4].position);
			}
		}
	
		/*!**********************************************************************
		 * \brief Given two sets of zoning data, calculate the distance in parameter space
		 * 
		 * \param i_new_rezone_data A pointer to an array of new rezone_union objects
		 * \param i_old_rezone_data A pointer to an array of old rezone_union objects
		 * 
		 * In order, the rezone_union objects must contain 1) a pointer to the element, 2) the total number of elements in the system, 3) min_size argument from rezone_minimize_ts, 4) max_size argument from rezone_minimize_ts, 5-n) positions of the zonal boundaries. This method is to be called by the GSL simulated annealing routine.
		 ************************************************************************/
		static double rezone_step_size (void *i_new_rezone_data, void *i_old_rezone_data) {
			datatype total = 0;
			bases::element <datatype> *element_ptr = ((rezone_union <datatype> *) i_new_rezone_data)->element_ptr;
			datatype new_positions [((rezone_union <datatype> *) i_new_rezone_data) [1].np + 1];
			for (int i = 0; i < ((rezone_union <datatype> *) i_new_rezone_data) [1].np + 1; ++i) {
				new_positions [i] = ((rezone_union <datatype> *) i_new_rezone_data) [i + 4].position;
			}
			datatype old_positions [((rezone_union <datatype> *) i_old_rezone_data) [1].np + 1];
			for (int i = 0; i < ((rezone_union <datatype> *) i_old_rezone_data) [1].np + 1; ++i) {
				old_positions [i] = ((rezone_union <datatype> *) i_old_rezone_data) [i + 4].position;
			}
			utils::messenger *messenger_ptr = element_ptr->messenger_ptr;
			for (int i = 1; i < messenger_ptr->get_np (); ++i) {
				total += (new_positions [i] - old_positions [i]) * (new_positions [i] - old_positions [i]);
			}
			return sqrt (total);
		}
	
		/*!**********************************************************************
		 * \brief Generate a new set of parameters from the old ones with a random number generator
		 * 
		 * \param r A pointer to a gnu scientific library random number generator
		 * \param i_rezone_data A pointer to an array of rezone_union objects
		 * \param step_size The datatype step_size in parameter space
		 * 
		 * In order, the rezone_union objects must contain 1) a pointer to the element, 2) the total number of elements in the system, 3) min_size argument from rezone_minimize_ts, 4) max_size argument from rezone_minimize_ts, 5-n) positions of the zonal boundaries. This method is to be called by the GSL simulated annealing routine.
		 ************************************************************************/
		static void rezone_generate_step (const gsl_rng *r, void *i_rezone_data, datatype step_size) {
			rezone_union <datatype> *rezone_data = (rezone_union <datatype> *) i_rezone_data;
			bases::element <datatype> *element_ptr = rezone_data->element_ptr;
			datatype positions [rezone_data [1].np + 1];
			for (int i = 0; i < rezone_data [1].np + 1; ++i) {
				positions [i] = rezone_data [i + 4].position;
			}
			utils::messenger *messenger_ptr = element_ptr->messenger_ptr;
			if (messenger_ptr->get_id () == 0) {
				// Generate a random radius that is less than or equal to the step size
				datatype radius = gsl_rng_uniform (r) * step_size;
				if (radius == 0.0) {
					messenger_ptr->skip_all ();
					return;
				}
				// Generate a random step for every zonal position and sum the total step size
				datatype xs [messenger_ptr->get_np () + 1];
				datatype total = 0.0;
				for (int i = 1; i < messenger_ptr->get_np (); ++i) {
					xs [i] = gsl_rng_uniform (r) * 2.0 - 1.0;
					total += xs [i] * xs [i];
				}
				// If possible, rescale the total step to be equal to the radius; otherwise, change nothing
				if (total == 0.0) {
					messenger_ptr->skip_all ();
					return;
				}
				total = sqrt (total);
				total /= radius;
				for (int i = 1; i < messenger_ptr->get_np (); ++i) {
					positions [i] = xs [i] / total + positions [i];
				}
				// Broadcast the new positions
				messenger_ptr->template bcast <datatype> (messenger_ptr->get_np () + 1, positions);
			} else {
				// Receive the new positions
				messenger_ptr->template bcast <datatype> (messenger_ptr->get_np () + 1, positions);
			}
			// Iterate forward through the data, checking that the mininum and maximum sizes are obeyed
			for (int i = 1; i < rezone_data [1].np; ++i) {
				// Check minimum size
				if (positions [i] < positions [i - 1] + rezone_data [2].position) {
					positions [i] = positions [i - 1] + rezone_data [2].position;
				}
				// Check maximum size
				if (positions [i] > positions [i - 1] + rezone_data [3].position) {
					positions [i] = positions [i - 1] + rezone_data [3].position;
				}
				rezone_data [i + 4].position = positions [i];
			}
			// Iterate backward through the positions, checking that the minimum and maximum sizes are obeyed
			for (int i = rezone_data [1].np - 1; i >= 1; --i) {
				// Check minimum size
				if (rezone_data [i + 4].position > rezone_data [i + 5].position - rezone_data [2].position) {
					rezone_data [i + 4].position = rezone_data [i + 5].position - rezone_data [2].position;
				}
				// Check maximum size
				if (rezone_data [i + 4].position < rezone_data [i + 5].position - rezone_data [3].position) {
					rezone_data [i + 4].position = rezone_data [i + 5].position - rezone_data [3].position;
				}
			}
			// Make sure we didn't accidentally run off the end of the position array
			if (rezone_data [4].position < rezone_data [5].position - rezone_data [3].position) {
				FATAL ("Unrealistic size constraints in rezone: max_size too small");
				throw 0;
			}
			if (rezone_data [4].position > rezone_data [5].position - rezone_data [2].position) {
				FATAL ("Unrealistic size constraints in rezone: min_size too large");
				throw 0;
			}
		}
	
		/*
			TODO Allow for upside-down elements
		*/
	};
} /* bases */

#endif /* end of include guard: ELEMENT_HPP_IUTSU4TQ */
