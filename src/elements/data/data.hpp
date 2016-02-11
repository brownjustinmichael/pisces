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
#include <ctime>
#include <chrono>
#include <unistd.h>

#include "linalg/utils.hpp"

#include "io/functors/slice.hpp"
#include "io/functors/product.hpp"
#include "io/input.hpp"
#include "io/formats/format.hpp"
#include "io/functors/functor.hpp"
#include "io/output.hpp"
#include "io/functors/average.hpp"
#include "io/formats/exceptions.hpp"
#include "io/formats/netcdf.hpp"
#include "plans/grids/grid.hpp"
#include "plans-transforms/implemented_transformer.hpp"

/*!**********************************************************************
 * \namespace data
 * 
 * \brief A namespace containing the data objects, which hold and output data
 ************************************************************************/
namespace data
{
	/*!**********************************************************************
	 * \brief A set of flags used when initializing data
	 ************************************************************************/
	enum initialize_element_flags {
		uniform_n = 0x01, //!< Copy the input array across the data in the n direction
		uniform_m = 0x02, //!< Copy the input array across the data in the m direction
		no_save = 0x04,
		no_variables = 0x08,
		profile = 0x80,
		vector = 0x10,
		vector2D = 0x10,
		vector3D = 0x20,
		corrector = 0x40
	};

	/*!**********************************************************************
	 * \brief A container object for the scalar data in the simulation
	 * 
	 * This object contains both the physical data and the output instructions.
	 ************************************************************************/
	class data
	{
	protected:
		int id; //!< The id of the current process
		io::parameters &params; //!< A reference to the parameters
		std::vector <grids::grid> grids; //!< A vector of the grid objects for each dimension

		static int mode; //!< The integer mode of the simulation (to prevent loading Chebyshev data into an even grid)
		
		std::vector <std::shared_ptr <io::output>> streams; //!< A vector of pointers to output streams
		std::vector <int> stream_conditions; //!< A vector of flags to match to output
		std::vector <bool> done; //!< A vector of booleans informing whether an output has been made this timestep
		
		std::shared_ptr <io::output> dump_stream; //!< A shared pointer to the dump output stream
		int dump_condition; //!< The flag condition of the data to be checked for dump
		bool dump_done; //!< A boolean stating whether the dump has happened
		
		std::vector <std::string> scalar_names; //!< The vector of variable names
		std::map <std::string, std::shared_ptr <grids::variable>> variables; //!< The string map of vectors containing the variable data
		double *weights; //!< Weights to be used for averaging over the data

		clock_t cbegin; //!< A beginning time, for timing
		clock_t cend; //!< An ending time, for timing 
		std::chrono::time_point <std::chrono::system_clock> tbegin; //!< A beginning time, for wall clock timing
		std::chrono::time_point <std::chrono::system_clock> tend; //!< An ending time, for wall clock timing

		double timing_cputime; //!< A double to keep track of cpu time
		double timing_walltime; //!< A double to keep track of wall time
		std::chrono::duration <double> timing_duration; //!< A duration object, for conversion from time_point to double
		
	public:
		int n_steps; //!< The number of steps the system has taken

		std::map <std::string, bool> is_corrector; //!< A map indicating whether an equation is a solver or a corrector (correctors take place only after all solvers are done)
		std::map <std::string, std::shared_ptr <plans::transforms::transformer>> transformers; //!< A map of shared pointers to the transformers associated with each variable

		double duration; //!< The total time elapsed in the simulation thus far
		double timestep; //!< The current timestep
		int flags; //!< A map of the simulation flags

		/**
		 * @param i_params A reference to the parameters object associated with the run
		 * @param i_id The integer identifier of this process, for output purposes
		 */
		data (io::parameters &i_params, int i_id = 0) : 
		id (i_id), 
		params (i_params) {
			flags = 0x0;
			duration = 0.0;
			cbegin=clock();
			tbegin=std::chrono::system_clock::now();
			timing_cputime = 0.0;
			timing_walltime = 0.0;
			timing_duration = std::chrono::duration <double>::zero ();
			n_steps = 0;
		}
		
		virtual ~data () {}

		/*!**********************************************************************
		 * \return The version of the class
		 ************************************************************************/
		static versions::version& version () {
			static versions::version version ("1.0.0.0");
			return version;
		}
		
		/*!**********************************************************************
		 * \brief Element iterator for iterating through the contained data
		 ************************************************************************/
		typedef std::vector <std::string>::iterator iterator;
	
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
		 * \param i_flags The flags indicating how to apply the initial conditions (uniform_n, uniform_m, or NULL)
		 * 
		 * This function initializes the dataset.
		 * 
		 * \return A pointer to the dataset
		 ************************************************************************/
		virtual grids::variable &initialize (std::string i_name, int i_flags = 0x00) {
			grids::variable &var = _initialize (i_name, i_flags);
			int output_flags = formats::all_d;
			if (i_flags & uniform_n) {
				output_flags = formats::m_profile;
			}

			if (dump_stream) {
				dump_stream->append <double> (i_name, var.ptr (), output_flags);
			}
			for (auto iter = streams.begin (); iter != streams.end (); ++iter)
			{
				(*iter)->append <double> (i_name, var.ptr (), output_flags);
			}
			scalar_names.push_back (i_name);
			std::sort (scalar_names.begin (), scalar_names.end ());
			is_corrector [i_name] = i_flags & corrector;
			return var;
		}
		
		/*!**********************************************************************
		 * \brief Check whether any streams can be output and output them
		 ************************************************************************/
		void output () {
			TRACE ("Output...");
			cend=clock();
			tend=std::chrono::system_clock::now();
			timing_cputime=((double) (cend - cbegin))/CLOCKS_PER_SEC;
			timing_duration=tend - tbegin;
			timing_walltime=timing_duration.count();
			for (int i = 0; i < (int) streams.size (); ++i) {
				streams [i]->to_file ();
			}
			if (dump_stream) dump_stream->to_file();
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
		
		/**
		 * @brief Return the string representation of the path to a file specified in a YAML Node
		 * @details The keys specified in the YAML Node can be "file" for the name of the file and "directory" for the directory in which the file should be located. If no directory is given, this will check the global parameters for "output.directory" to use instead. The output file can have two format-specified integers as %i and %%i. The first will be replaced with the id number of the process. The second will be replaced with the contents of "output.number" or zero.
		 * 
		 * @param i_params The Yaml Node from which the file name is to be constructed
		 * @return The string representation of the path to a file specified in the given YAML Node
		 */
		std::string file_from (YAML::Node i_params) {
			if (!(i_params.IsDefined ())) return "";

			std::string file_format = "";
			if (i_params ["file"].IsDefined ()) file_format = i_params ["file"].as <std::string> ();

			if (file_format == "") return "";

			if (i_params ["directory"].IsDefined ()) {
				file_format = i_params ["directory"].as <std::string> () + file_format;
			} else if (params ["output.directory"].IsDefined ()) {
				file_format = params ["output.directory"].as <std::string> () + file_format;
			}

			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), id);

			file_format = buffer;

			snprintf (buffer, file_format.size () * 2, file_format.c_str (), params ["output.number"].as <int> ());

			return buffer;
		}

		/**
		 * @brief Set up the data object from an input stream
		 * @details Read in the input stream into the data object. Warnings will be raised if any variables appear to be missing, but anything absent will be assumed to be zero, and the computation will commence
		 * 
		 * @param input_stream The input stream to read in
		 * @param state The state of the data into which the data should be read in
		 * 
		 * @return A pointer to the input stream parameter
		 */
		virtual std::shared_ptr <io::input> setup (std::shared_ptr <io::input> input_stream, int state = real_real) {
			for (data::iterator iter = begin (); iter != end (); ++iter) {
				input_stream->append (*iter, (*this) (*iter, state));
				(*this) [*iter].last_update = state;
			}
			
			input_stream->append ("t", &duration, formats::scalar);
			input_stream->append ("dt", &timestep, formats::scalar);
			input_stream->append ("step", &n_steps, formats::scalar);

			try {
				// Read from the input into the element variables
				input_stream->from_file ();

			} catch (formats::exceptions::bad_variables &except) {
				// Send a warning if there are missing variables in the input
				WARN (except.what ());
			}

			for (data::iterator iter = begin (); iter != end (); ++iter) {
				if (transformers [*iter]) transformers [*iter]->update ();
			}

			return input_stream;
		}

		/**
		 * @brief Set up the data from a YAML Node of the input parameters and a grid
		 * @details This method is added for a more convenient call to the setup tool.
		 * 
		 * @param input_params The YAML Node describing the parameters definining the input file
		 * @param grid a data_grid object that describes the dimensions of the input file
		 * 
		 * @return A shared pointer to the input object
		 */
		template <class format>
		std::shared_ptr <io::input> setup_from (YAML::Node input_params, const formats::data_grid &grid) {
			if (file_from (input_params) == "") return std::shared_ptr <io::input> ();

			std::shared_ptr <io::input> input_stream (new io::formatted_input <format> (grid, file_from (input_params)));

			return setup (input_stream);
		}

		/**
		 * @brief Add a global attribute to each stream associated with the data object
		 * @details Outputs can have string attributes associated with them to allow for metadata storage.
		 * 
		 * @param name The key of the associated attribute
		 * @param attribute The string attribute to add to the output
		 */
		void add_stream_attribute (std::string name, std::string attribute) {
			for (auto iter = streams.begin (); iter < streams.end (); ++iter)
			{
				(*iter)->add_global_attribute (name, attribute);
			}
		}
	
		/*!**********************************************************************
		 * \brief Given an output_stream, append all relevant outputs
		 * 
		 * \param output A shared_ptr to the output object
		 * @param state The integer state from which the output should be read
		 * \param flags The integer flags to specify the type of output (transform stream, normal_stream)
		 * 
		 * Given an output object, prepare it to output all the scalar fields tracked directly by the element. Make it one of the output streams called during output. If flags is normal_stream, output the variables when in Cartesian space, else if flags is transform_stream, output with horizontal grid in Fourier space, but vertical grid in Cartesian space.
		 ************************************************************************/
		std::shared_ptr <io::output> setup_output (std::shared_ptr <io::output> output, int state = real_real, int flags = 0x00) {
			// Iterate through the scalar fields and append them to the variables which the output will write to file
			TRACE ("Setting up output stream...");

			if (!(params ["output.output"].as <bool> ())) return std::shared_ptr <io::output> ();

			if (!(flags & no_variables)) {
				if (flags & profile) {
					for (data::iterator iter = begin (); iter != end (); ++iter) {
						TRACE ("Appending " << *iter << " profile to output...");
						output->append<double>(*iter, this->output_prof(*iter), formats::m_profile);
						output->append<double>(*iter + "_deriv", this->output_deriv_prof(*iter), formats::m_profile);
					}
				} else {
					for (data::iterator iter = begin (); iter != end (); ++iter) {
						TRACE ("Appending " << *iter << " to output...");
						output->append (*iter, (*this) (*iter, state));
					}
				}
			}
			
			output->append ("t", &duration, formats::scalar);
			output->append ("dt", &timestep, formats::scalar);
			output->append ("cputime", &timing_cputime, formats::scalar);
			output->append ("walltime", &timing_walltime, formats::scalar);
			output->append("step", &n_steps, formats::scalar);

			output->add_global_attribute ("params", params.string ());
			output->add_global_attribute ("data_version", version ());
			output->add_global_attribute ("output_version", output->version ());

			char buff [1000];
			gethostname (buff, sizeof(buff));
			std::string host = buff;
			output->add_global_attribute("hostname", host);

			getcwd (buff, sizeof(buff));
			host = buff;		
			output->add_global_attribute("cwd", host);

			// Check the desired output time and save the output object in the appropriate variable
			streams.push_back (output);

			return output;
		}

		/**
		 * @brief Set up an output stream from a YAML Node
		 * @details Given a YAML Node, which can have various keys, construct an output stream. These keys can be "file" for the file name, "timed" for whether to output every set amount of time or number of time steps, "timed_every" for the duration of time between outputs if "timed" was true, and "every" for the number of time steps between outputs if "timed" was false.
		 * 
		 * @param output_params The YAML Node associated with the output stream
		 * @param grid The data_grid containing the information about the desired dimensions to output
		 * @param state The state to output
		 * @param flags Any flags associated with the output
		 * @return A shared pointer to the newly constructed output stream
		 */
		template <class format>
		std::shared_ptr <io::output> setup_output_from (YAML::Node output_params, const formats::data_grid &grid, int state = real_real, int flags = 0x00) {
			// Iterate through the scalar fields and append them to the variables for which the input will search
			TRACE ("Setting up output from YAML Node...");
			if (file_from (output_params) == "") return std::shared_ptr <io::output> ();

			std::shared_ptr <io::output> output_stream;
			if (output_params ["timed"].IsDefined () && output_params ["timed"].as <bool> ()) {
				output_stream.reset (new io::timed_appender_output <double, formats::netcdf> (grid, file_from (output_params), duration, output_params ["timed_every"].as <double> ()));
			} else {
				output_stream.reset (new io::appender_output <formats::netcdf> (grid, file_from (output_params), &duration, output_params ["every"].as <int> ()));
			}
			return setup_output (output_stream, state, flags);
		}
		
		/*!**********************************************************************
		 * \brief Set up the dump stream
		 * 
		 * \param output_ptr The shared pointer to the output stream to use as the dump
		 * \param flags The conditions to be met for dump (e.g. transformed_vertical | transformed_horizontal)
		 ************************************************************************/
		virtual void setup_dump (std::shared_ptr <io::output> output_ptr, int flags = 0x00) {
			// Iterate through the scalar fields and append them to the variables which the output will write to file
			output_ptr->append <double> ("t", &duration, formats::scalar);
			output_ptr->append <double> ("dt", &timestep, formats::scalar);
			output_ptr->append <const int> ("mode", &(get_mode ()), formats::scalar);
			output_ptr->append("step", &n_steps, formats::scalar);

			// Check the desired output time and save the output object in the appropriate variable
			dump_stream = output_ptr;
			dump_done = false;
			dump_condition = flags;
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
		
		/*!**********************************************************************
		 * \brief Given an integer index to one of the scalar variables, return a pointer to that dataset
		 * 
		 * \param name The integer name to index
		 * 
		 * \return A pointer to the indexed data
		 ************************************************************************/
		grids::variable &operator[] (const std::string &name) {
			return *variables [name];
		}
		
		/*!**********************************************************************
		 * \brief Given an integer name to one of the scalar variables and an index within that dataset, return a pointer to that dataset at the given index
		 * 
		 * \param name The integer name to index
		 * @param state The integer state of the variable from which the data should be directed
		 * \param index The integer index inside the named variable
		 * 
		 * \return A pointer to the indexed data
		 ************************************************************************/
		double *operator() (const std::string &name, int state = 0, int index = 0) {
			if (variables.find (name) != variables.end () && (int) variables.size () != 0) {
				return variables [name]->ptr (state) + index * variables [name]->dims ();
			}
			ERROR ("Variable " << name << " is undefined");
			throw 500;
		}

		/**
		 * @brief Construct a functor that returns the profile of a given string variable
		 * 
		 * @param variable The string variable from which the profile should be calculated
		 * @return The profile functor, which should be added to an output stream
		 */
		virtual std::shared_ptr <functors::functor> output_prof (std::string variable) = 0;

		/**
		 * @brief Construct a functor that returns the profile of the derivative a given string variable
		 * 
		 * @param variable The string variable from which the profile should be calculated
		 * @return The profile functor, which should be added to an output stream
		 */
		virtual std::shared_ptr <functors::functor> output_deriv_prof (std::string variable) = 0;

	protected:
		/*!**********************************************************************
		 * \brief Initialize a new dataset in the data object
		 ************************************************************************/
		virtual grids::variable &_initialize (std::string i_name, int i_flags = 0x00) = 0;
	};
} /* data */

#endif /* end of include guard: DATA_HPP_574F2B9E */
