/*!**********************************************************************
 * \file data.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-10-09.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DATA_HPP_574F2B9E
#define DATA_HPP_574F2B9E

#include <string>
#include <map>
#include <memory>
#include <algorithm>

#include "linalg/utils.hpp"

#include "io/input.hpp"
#include "io/formats/format.hpp"
#include "io/functors/functor.hpp"
#include "io/output.hpp"
#include "io/functors/average.hpp"
#include "io/formats/exceptions.hpp"
#include "io/formats/netcdf.hpp"
#include "plans/grids/grid.hpp"
#include "plans-transforms/implemented_transformer.hpp"

namespace data
{
	/*!**********************************************************************
	 * \brief A set of flags used when initializing data
	 ************************************************************************/
	enum initialize_element_flags {
		uniform_n = 0x01, //!< Copy the input array across the data in the n direction
		uniform_m = 0x02, //!< Copy the input array across the data in the m direction
		no_save = 0x04,
	};
	
	/*!**********************************************************************
	 * \brief A container object for the scalar data in the simulation
	 * 
	 * This object contains both the physical data and the output instructions.
	 ************************************************************************/
	template <class datatype>
	class data
	{
	protected:
		static int mode; //!< The integer mode of the simulation (to prevent loading Chebyshev data into an even grid)
		
		std::vector <std::shared_ptr <io::output>> streams; //!< A vector of pointers to output streams
		std::vector <int> stream_conditions; //!< A vector of flags to match to output
		std::vector <bool> done; //!< A vector of booleans informing whether an output has been made this timestep
		
		std::shared_ptr <io::output> dump_stream;
		int dump_condition;
		bool dump_done;
		
	public:
		datatype duration; //!< The total time elapsed in the simulation thus far
		datatype timestep; //!< The current timestep
		std::map <std::string, int> flags; //!< A map of the simulation flags
		data () {
			flags ["state"] = 0x0;
			duration = 0.0;
		}
		
		virtual ~data () {}
		
		/*!**********************************************************************
		 * \brief Given an integer index to one of the scalar variables, return a pointer to that dataset
		 * 
		 * \param name The integer name to index
		 * 
		 * \return A pointer to the indexed data
		 ************************************************************************/
		datatype *operator[] (const std::string &name) {
			return &(scalars [name] [0]);
		}
		
		/*!**********************************************************************
		 * \brief Element iterator for iterating through the contained data
		 ************************************************************************/
		typedef typename std::vector <std::string>::iterator iterator;
	
		/*!**********************************************************************
		 * \brief Generate an iterator to iterate through the data
		 * 
		 * \return Beginning iterator
		 ************************************************************************/
		iterator begin () {
			return scalar_names.begin ();
		}
	
		/*!**********************************************************************
		 * \brief Generate a finishing iterator for comparison
		 * 
		 * \return Ending iterator
		 ************************************************************************/
		iterator end () {
			return scalar_names.end ();
		}
		
		/*!**********************************************************************
		 * \brief Given an integer name to one of the scalar variables and an index within that dataset, return a pointer to that dataset at the given index
		 * 
		 * \param name The integer name to index
		 * \param index The integer index inside the named variable
		 * 
		 * \return A pointer to the indexed data
		 ************************************************************************/
		datatype *operator() (const std::string &name, const int index = 0) {
			if (scalars.find (name) != scalars.end () && (int) scalars [name].size () != 0) {
				return &(scalars [name] [index]);
			}
			return NULL;
		}
		
		/*!**********************************************************************
		 * \brief Get the mode of the simulation
		 * 
		 * \return A reference to the integer mode of the simulation
		 ************************************************************************/
		const int &get_mode () {
			return mode;
		}
		
		/*!**********************************************************************
		 * \brief Initialize a dataset for the simulation
		 * 
		 * \param i_name The integer name to give the new dataset
		 * \param i_str The string name of the new dataset
		 * \param initial_conditions The initial values of the dataset
		 * \param i_flags The flags indicating how to apply the initial conditions (uniform_n, uniform_m, or NULL)
		 * 
		 * This function initializes the dataset.
		 * 
		 * \return A pointer to the dataset
		 ************************************************************************/
		virtual datatype *initialize (std::string i_name, datatype* initial_conditions = NULL, int i_flags = 0x00) {
			flags [i_name] = 0x00;
			datatype *ptr = _initialize (i_name, initial_conditions, i_flags);
			if (dump_stream) {
				dump_stream->template append <datatype> (i_name, ptr);
			}
			scalar_names.push_back (i_name);
			std::sort (scalar_names.begin (), scalar_names.end ());
			return ptr;
		}
		
		/*!**********************************************************************
		 * \brief Check whether any streams can be output and output them
		 ************************************************************************/
		void check_streams (int i_flags = 0x00) {
			TRACE ("Transforming...");
			for (int i = 0; i < (int) streams.size (); ++i) {
				if (done [i]) continue;
				if (stream_conditions [i] & transformed_vertical) {
					if (!(i_flags & transformed_vertical)) {
						continue;
					}
				} else {
					if ((i_flags & transformed_vertical)) {
						continue;
					}
				}
				if (stream_conditions [i] & transformed_horizontal) {
					if (!(i_flags & transformed_horizontal)) {
						continue;
					}
				} else {
					if ((i_flags & transformed_horizontal)) {
						continue;
					}
				}
				streams [i]->to_file ();
				done [i] = true;
			}
			if (dump_stream && !dump_done) {
				if (dump_condition & transformed_vertical) {
					if (!(i_flags & transformed_vertical)) {
						return;
					}
				} else {
					if ((i_flags & transformed_vertical)) {
						return;
					}
				}
				if (dump_condition & transformed_horizontal) {
					if (!(i_flags & transformed_horizontal)) {
						return;
					}
				} else {
					if ((i_flags & transformed_horizontal)) {
						return;
					}
				}
				dump_stream->to_file ();
				dump_done = true;
			}
		}
		
		/*!**********************************************************************
		 * \brief Every timestep, reset the output stream conditions
		 ************************************************************************/
		void reset () {
			for (int i = 0; i < (int) done.size (); ++i) {
				done [i] = false;
			}
			dump_done = false;
		}
		
		/*!**********************************************************************
		 * \brief Given an input stream, load all the relevant data into the data object
		 * 
		 * \param input_stream_ptr A pointer to an input object
		 ************************************************************************/
		void setup (io::input *input_ptr) {
			// Iterate through the scalar fields and append them to the variables for which the input will search
			for (data::iterator iter = begin (); iter != end (); ++iter) {
				input_ptr->template append <datatype> (*iter, (*this) (*iter));
			}
			
			input_ptr->template append <datatype> ("t", &duration, io::scalar);
			int mode;
			input_ptr->template append <int> ("mode", &mode, io::scalar);
			
			try {
				// Read from the input into the element variables
				input_ptr->from_file ();
			} catch (io::formats::exceptions::bad_variables &except) {
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
		virtual void setup_output (std::shared_ptr <io::output> output_ptr, int flags = 0x00) {
			// Iterate through the scalar fields and append them to the variables which the output will write to file
			for (data::iterator iter = begin (); iter != end (); ++iter) {
				output_ptr->template append <datatype> (*iter, (*this) (*iter));
			}
			
			output_ptr->template append <datatype> ("t", &duration, io::scalar);
			output_ptr->template append <const int> ("mode", &(get_mode ()), io::scalar);

			// Check the desired output time and save the output object in the appropriate variable
			if (!(flags & no_save)) {
				streams.push_back (output_ptr);
				done.push_back (false);
				stream_conditions.push_back (flags);
			}
		}
		
		virtual void setup_dump (std::shared_ptr <io::output> output_ptr, int flags = 0x00) {
			// Iterate through the scalar fields and append them to the variables which the output will write to file
			output_ptr->template append <datatype> ("t", &duration, io::scalar);
			output_ptr->template append <const int> ("mode", &(get_mode ()), io::scalar);

			// Check the desired output time and save the output object in the appropriate variable
			dump_stream = output_ptr;
			dump_done = false;
			dump_condition = flags;
		}
		
		/*!**********************************************************************
		 * \brief Given an output stream, prepare the stat output
		 * 
		 * \param output_ptr A shared_ptr to the output object
		 * \param flags The integer flags to specify the time of output
		 ************************************************************************/
		virtual void setup_stat (std::shared_ptr <io::output> output_ptr, int flags = 0x00) {
			// Also prepare to output the total simulated time and geometry mode
			output_ptr->template append <datatype> ("t", &duration, io::scalar);
			output_ptr->template append <datatype> ("dt", &timestep, io::scalar);
			// Check the desired output time and save the output object in the appropriate variable
			if (!(flags & no_save)){
				streams.push_back (output_ptr);
				stream_conditions.push_back (0x00);
				done.push_back (false);
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
			// normal_profiles.push_back (output_ptr);
		}
	
		/*
			TODO These three methods are related in purpose but varied in design; they should be unified
		*/
		
	protected:
		std::vector <std::string> scalar_names;
		std::map <std::string, std::vector <datatype>> scalars;
		
		/*!**********************************************************************
		 * \brief Initialize a new dataset in the data object
		 ************************************************************************/
		virtual datatype *_initialize (std::string i_name, datatype* initial_conditions = NULL, int i_flags = 0x00) = 0;
	};
	
	template <class datatype>
	class implemented_data : public data <datatype>
	{
	protected:
		using data <datatype>::iterator;
		std::shared_ptr <grids::grid <datatype>> grid_n, grid_m;
		int n, m;
		
	public:
		using data <datatype>::duration;
		implemented_data (grids::axis *i_axis_n, grids::axis *i_axis_m, int name = 0, std::string dump_file = "", std::string dump_directory = "./", int dump_every = 1) : grid_n (std::shared_ptr <grids::grid <datatype>> (new typename grids::horizontal::grid <datatype> (i_axis_n))), grid_m (std::shared_ptr <grids::grid <datatype>> (new typename grids::vertical::grid <datatype> (i_axis_m))), n (grid_n->get_n ()), m (grid_m->get_n ()) {
			// Set up output
			const io::data_grid o_grid = io::data_grid::two_d (n, m, 0, 0, 0, 0);
			
			std::shared_ptr <io::output> dump_stream;
			if (dump_file != "") {
				std::string file_format = dump_directory + dump_file;
				char buffer [file_format.size () * 2];
				snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);

				dump_stream.reset (new io::replace_output <io::formats::two_d::netcdf> (o_grid, buffer, dump_every));
				this->setup_dump (dump_stream);
			}
			
			this->initialize ("x");
			this->initialize ("z");
			
			/*
				TODO Clean dump file generation
			*/
		}
		
		virtual ~implemented_data () {}
		
		datatype *operator () (const std::string name, const int i = 0, const int j = 0) {
			return data <datatype>::operator() (name, i * grid_m->get_ld () + j); 
		}
		
		virtual void setup_profile (std::shared_ptr <io::output> output_ptr, int flags = 0x00) {
			typedef typename std::map <std::string, std::vector <datatype> >::iterator iterator;
			for (iterator iter = data <datatype>::scalars.begin (); iter != data <datatype>::scalars.end (); ++iter) {
				output_ptr->template append <datatype> (iter->first, (*this) (iter->first));
				output_ptr->template append <datatype> ("rms_" + iter->first, std::shared_ptr <io::functors::functor> (new io::functors::root_mean_square_functor <datatype> ((*this) (iter->first), n, m)));
				output_ptr->template append <datatype> ("avg_" + iter->first, std::shared_ptr <io::functors::functor> (new io::functors::average_functor <datatype> ((*this) (iter->first), n, m)));
			}
			output_ptr->template append <datatype> ("t", &(duration), io::scalar);
			output_ptr->template append <const int> ("mode", &(data <datatype>::get_mode ()), io::scalar);

			data <datatype>::setup_profile (output_ptr, flags);
		}
	
	protected:
		virtual datatype *_initialize (std::string name, datatype* initial_conditions = NULL, int i_flags = 0x00) {
			TRACE ("Initializing " << name << "...");
			// Size allowing for real FFT buffer
			data <datatype>::scalars [name].resize (grid_n->get_ld () * m, 0.0);
			if (name == "x") {
				for (int j = 0; j < m; ++j) {
					linalg::copy (n, &((*grid_n) [0]), (*this) (name, 0, j), 1, m);
				}
			} else if (name == "z") {
				for (int i = 0; i < n; ++i) {
					linalg::copy (m, &((*grid_m) [0]), (*this) (name, i));
				}
			} else {
				if (initial_conditions) {
					if ((i_flags & uniform_m) && (i_flags & uniform_n)) {
						for (int i = 0; i < n; ++i) {
							for (int j = 0; j < m; ++j) {
								*(*this) (name, i, j) = *initial_conditions;
							}
						}
					} else if (i_flags & uniform_m) {
						for (int j = 0; j < m; ++j) {
							linalg::copy (n, initial_conditions, (*this) (name, 0, j), 1, m);
						}
					} else if (i_flags & uniform_n) {
						for (int i = 0; i < n; ++i) {
							linalg::copy (m, initial_conditions, (*this) (name, i));
						}
					} else {
						linalg::copy (n * m, initial_conditions, (*this) (name));
					}
				}
			}
			TRACE ("Done.");
			return &(data <datatype>::scalars [name] [0]);
		}
	};
} /* data */

#endif /* end of include guard: DATA_HPP_574F2B9E */
