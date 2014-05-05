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
#include <gsl/gsl_siman.h>

namespace bases
{
	template <class datatype>
	class element;
	
	template <class datatype>
	union rezone_union {
		bases::element <datatype> *element_ptr;
		int np;
		datatype position;
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
		virtual datatype *initialize (int i_name, std::string i_str, datatype* initial_conditions = NULL, int i_flags = 0x00) {
			element_flags [i_name] = 0x00;
			for (int i = 0; i < dimensions; ++i) {
				if (!grids [i]) {
					grids [i] = generate_grid (&axes [i], i);
				}
			}
			scalar_names [i_name] = i_str;
			return _initialize (i_name, initial_conditions, i_flags);
		}
		
		virtual std::shared_ptr <grid <datatype>> generate_grid (bases::axis *axis, int index = -1) = 0;
		
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
		
		virtual const int &get_mode () = 0;
	
		void setup (io::input *input_stream) {
			// Set up input
			typedef typename std::map <int, std::vector <datatype> >::iterator iterator;
			for (iterator iter = scalars.begin (); iter != scalars.end (); ++iter) {
				input_stream->template append <datatype> (scalar_names [iter->first], ptr (iter->first));
			}
			input_stream->template append_scalar <datatype> ("t", &duration);
			int mode;
			input_stream->template append_scalar <int> ("mode", &mode);

			try {
				input_stream->from_file ();
			} catch (exceptions::io::bad_variables &except) {
				WARN (except.what ());
			}
			

			if (mode != get_mode ()) {
				FATAL ("Loading simulation in different mode: " << mode << " instead of " << get_mode ());
				throw 0;
			}
		}
		
		void setup_output (std::shared_ptr <io::output> output_stream, int flags = 0x00) {
			/*
				TODO This breaks if a pointer to an existing stream is passed, rather than a new instance
			*/
			typedef typename std::map <int, std::vector <datatype> >::iterator iterator;
			for (iterator iter = scalars.begin (); iter != scalars.end (); ++iter) {
				output_stream->template append <datatype> (scalar_names [iter->first], ptr (iter->first));
			}
			output_stream->template append_scalar <datatype> ("t", &duration);
			output_stream->template append_scalar <const int> ("mode", &(get_mode ()));

			if (flags & transform_output) {
				transform_stream = output_stream;
			} else if (flags & normal_output) {
				normal_stream = output_stream;
			}
		}
		
		virtual void setup_profile (std::shared_ptr <io::output> output_stream, int flags = 0x00) {
			normal_profiles.push_back (output_stream);
		}
		
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
		
		virtual std::shared_ptr <io::virtual_dump> make_dump (int flags = 0x00) = 0;
		
		virtual std::shared_ptr <io::virtual_dump> make_rezoned_dump (datatype *positions, io::virtual_dump *dump, int flags = 0x00) = 0;
		
		virtual void get_zoning_positions (datatype *positions) = 0;
		
		/*!**********************************************************************
		 * \param n_tries the number of attempts before stepping
		 * \param iters_fixed_t the integer number of iterations at each temperature
		 * \param step_size the real maximum step size
		 * \param k a real boltzmann constant
		 * \param t_initial the real initial temperature
		 * \param mu_t the real damping factor for temperature
		 * \param t_min the real damping factor parameter for temperature
		 ************************************************************************/
		virtual std::shared_ptr <io::virtual_dump> rezone_minimize_ts (datatype * positions, datatype min_size, datatype max_size, int n_tries = 20, int iters_fixed_t = 1000, datatype step_size = 1.0, datatype k = 1.0, datatype t_initial = 0.008, datatype mu_t = 1.003, datatype t_min = 2.0e-6) {
			TRACE ("Rezoning...");
			write_transform_data ();
			transform (inverse_horizontal | inverse_vertical);
			read_transform_data ();

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
			
			DEBUG ("Old timestep: " << timestep << " New timestep: " << -rezone_calculate_ts (rezone_data));
			if (rezone_calculate_ts (rezone_data) < -timestep) {
				for (int i = 0; i < messenger_ptr->get_np () + 1; ++i) {
					positions [i] = rezone_data [i + 4].position;
				}
			}
			
			DEBUG ("New positions: " << positions [messenger_ptr->get_id ()] << " " << positions [messenger_ptr->get_id () + 1]);
			
			return make_rezoned_dump (positions, &*(make_dump ()));
		}
		
		static double rezone_calculate_ts (void *i_rezone_data) {
			bases::element <datatype> *element_ptr = ((rezone_union <datatype> *) i_rezone_data)->element_ptr;
			datatype positions [((rezone_union <datatype> *) i_rezone_data) [1].np + 1];
			bases::messenger *messenger_ptr = element_ptr->messenger_ptr;
			for (int i = 0; i < ((rezone_union <datatype> *) i_rezone_data) [1].np + 1; ++i) {
				positions [i] = ((rezone_union <datatype> *) i_rezone_data) [i + 4].position;
			}
			double timestep = element_ptr->calculate_min_timestep (&*(element_ptr->make_rezoned_dump (positions, &*(element_ptr->rezone_dump), profile_only)));
			messenger_ptr->min (&timestep);
			return -timestep;
		}
		
		static void print_rezone (void *i_rezone_data) {
			for (int i = 0; i < ((rezone_union <datatype> *) i_rezone_data) [1].np + 1; ++i) {
				printf (" %f ", ((rezone_union <datatype> *) i_rezone_data) [i + 4].position);
			}
		}
		
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
			bases::messenger *messenger_ptr = element_ptr->messenger_ptr;
			for (int i = 1; i < messenger_ptr->get_np (); ++i) {
				total += (new_positions [i] - old_positions [i]) * (new_positions [i] - old_positions [i]);
			}
			return sqrt (total);
		}
		
		static void rezone_generate_step (const gsl_rng *r, void *i_rezone_data, double step_size) {
			rezone_union <datatype> *rezone_data = (rezone_union <datatype> *) i_rezone_data;
			bases::element <datatype> *element_ptr = rezone_data->element_ptr;
			datatype positions [rezone_data [1].np + 1];
			for (int i = 0; i < rezone_data [1].np + 1; ++i) {
				positions [i] = rezone_data [i + 4].position;
			}
			bases::messenger *messenger_ptr = element_ptr->messenger_ptr;
			if (messenger_ptr->get_id () == 0) {
				datatype radius = gsl_rng_uniform (r) * step_size;
				if (radius == 0.0) {
					int flags = mpi_skip;
					messenger_ptr->check_all (&flags);
					return;
				}
				datatype xs [messenger_ptr->get_np () + 1];
				datatype total = 0.0;
				for (int i = 1; i < messenger_ptr->get_np (); ++i) {
					xs [i] = gsl_rng_uniform (r) * 2.0 - 1.0;
					total += xs [i] * xs [i];
				}
				if (total == 0.0) {
					int flags = mpi_skip;
					messenger_ptr->check_all (&flags);
					return;
				}
				total = sqrt (total);
				total /= radius;
				for (int i = 1; i < messenger_ptr->get_np (); ++i) {
					positions [i] = xs [i] / total + positions [i];
				}
				messenger_ptr->template bcast <datatype> (messenger_ptr->get_np () + 1, positions);
			} else {
				messenger_ptr->template bcast <datatype> (messenger_ptr->get_np () + 1, positions);
			}
			for (int i = 1; i < rezone_data [1].np; ++i) {
				if (positions [i] < positions [i - 1] + rezone_data [2].position) {
					positions [i] = positions [i - 1] + rezone_data [2].position;
				}
				if (positions [i] > positions [i - 1] + rezone_data [3].position) {
					positions [i] = positions [i - 1] + rezone_data [3].position;
				}
				rezone_data [i + 4].position = positions [i];
			}
			for (int i = rezone_data [1].np - 1; i >= 1; --i) {
				if (rezone_data [i + 4].position > rezone_data [i + 5].position - rezone_data [2].position) {
					rezone_data [i + 4].position = rezone_data [i + 5].position - rezone_data [2].position;
				}
				if (rezone_data [i + 4].position < rezone_data [i + 5].position - rezone_data [3].position) {
					rezone_data [i + 4].position = rezone_data [i + 5].position - rezone_data [3].position;
				}
			}
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
		
		/*!**********************************************************************
		 * \brief The main function call of the class
		 * 
		 * This method tells the element to begin the main run of the simulation.
		 * It runs through all the specified plans in the appropriate order, and 
		 * updates the values as necessary. Output, if desired, is specified by 
		 * the output streams.
		 ************************************************************************/
		virtual void run (int &n_steps, int max_steps);
		
		/*!*******************************************************************
		 * \brief Output all information to a dump file in the current directory
		 * 
		 * This failsafe method is designed to create a snapshot immediately
		 * before the simulation crashed.
		 *********************************************************************/
		inline virtual void failsafe () {
			failsafe_dump->to_file ();
		}
		
		std::vector <bases::axis> axes;
		std::vector <std::shared_ptr <grid <datatype>>> grids; //!< A vector of shared pointers to the collocation grids
		std::shared_ptr <io::virtual_dump> rezone_dump;
		messenger* messenger_ptr; //!< A pointer to the messenger object
	protected:
		int name, dimensions; //!< An integer representation of the element, to be used in file output
		io::parameters& params; //!< The map that contains the input parameters
		
		datatype duration; //!< The datatype total simulated time
		datatype timestep; //!< The datatype timestep length

		std::map <int, std::vector <datatype> > scalars; //!< A map of scalar vectors
		std::map <int, std::string> scalar_names;
		std::map <int, int> element_flags; //!< A map of integer flags
		
		std::shared_ptr <io::output> failsafe_dump; //!< An implementation to dump in case of failure
		std::shared_ptr <io::output> normal_stream; //!< An implementation to output in normal space
		std::shared_ptr <io::output> transform_stream; //!< An implementation to output in transform space
		std::vector <std::shared_ptr <io::output>> normal_profiles; //!< An implementation to output in transform space
		
		
		std::vector <int> solver_keys; //!< A vector of integer keys to the solvers map
		std::map<int, std::shared_ptr<solver <datatype>>> solvers; //!< A vector of shared pointers to the matrix solvers

	private:
		std::vector <int> transforms; //!< A vector of integer keys to the transform maps
		std::map <int, std::shared_ptr <master_transform <datatype>>> master_transforms;
	};
} /* bases */

#endif /* end of include guard: ELEMENT_HPP_IUTSU4TQ */
