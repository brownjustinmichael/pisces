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
	template <class datatype>
	class data
	{
	protected:
		int id;
		io::parameters &params;
		std::vector <grids::grid <datatype>> grids;

		static int mode; //!< The integer mode of the simulation (to prevent loading Chebyshev data into an even grid)
		
		std::vector <std::shared_ptr <io::output>> streams; //!< A vector of pointers to output streams
		std::vector <int> stream_conditions; //!< A vector of flags to match to output
		std::vector <bool> done; //!< A vector of booleans informing whether an output has been made this timestep
		
		std::shared_ptr <io::output> dump_stream; //!< A shared pointer to the dump output stream
		int dump_condition; //!< The flag condition of the data to be checked for dump
		bool dump_done; //!< A boolean stating whether the dump has happened
		
		std::vector <std::string> scalar_names; //!< The vector of variable names
		std::map <std::string, std::shared_ptr <grids::variable <datatype>>> variables; //!< The string map of vectors containing the variable data
		datatype *weights;

		clock_t cbegin, cend;
		std::chrono::time_point <std::chrono::system_clock> tbegin, tend;

		double timing_cputime = 0.0;
		double timing_walltime = 0.0;
		std::chrono::duration <double> timing_duration = std::chrono::duration <double>::zero ();
		
	public:
		int n_steps = 0;

		std::map <std::string, bool> is_corrector;
		std::map <std::string, std::shared_ptr <plans::transforms::transformer <datatype>>> transformers;

		datatype duration; //!< The total time elapsed in the simulation thus far
		datatype timestep; //!< The current timestep
		int flags; //!< A map of the simulation flags
		data (io::parameters &i_params, int i_id = 0) : id (i_id), params (i_params) {
			flags = 0x0;
			duration = 0.0;
			cbegin=clock();
			tbegin=std::chrono::system_clock::now();
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
		 * \brief Given an integer index to one of the scalar variables, return a pointer to that dataset
		 * 
		 * \param name The integer name to index
		 * 
		 * \return A pointer to the indexed data
		 ************************************************************************/
		grids::variable <datatype> &operator[] (const std::string &name) {
			return *variables [name];
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
		datatype *operator() (const std::string &name, int state = 0, int index = 0) {
			if (variables.find (name) != variables.end () && (int) variables.size () != 0) {
				return variables [name]->ptr (state) + index * variables [name]->dims ();
			}
			ERROR ("Variable " << name << " is undefined");
			throw 500;
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
		 * \param initial_conditions The initial values of the dataset
		 * \param i_flags The flags indicating how to apply the initial conditions (uniform_n, uniform_m, or NULL)
		 * 
		 * This function initializes the dataset.
		 * 
		 * \return A pointer to the dataset
		 ************************************************************************/
		virtual grids::variable <datatype> &initialize (std::string i_name, int i_flags = 0x00) {
			grids::variable <datatype> &var = _initialize (i_name, i_flags);
			if (dump_stream) {
				dump_stream->template append <datatype> (i_name, var.ptr ());
			}
			for (auto iter = streams.begin (); iter != streams.end (); ++iter)
			{
				(*iter)->template append <datatype> (i_name, var.ptr ());
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
			dump_stream->to_file();
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

		template <class format>
		std::shared_ptr <io::input> setup_from (YAML::Node input_params, const formats::data_grid &grid) {
			if (file_from (input_params) == "") return NULL;

			std::shared_ptr <io::input> input_stream (new io::formatted_input <format> (grid, file_from (input_params)));

			return setup (input_stream);
		}

		void add_stream_attribute (std::string name, std::string attribute) {
			for (auto iter = streams.begin (); iter < streams.end (); ++iter)
			{
				(*iter)->add_global_attribute (name, attribute);
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
		std::shared_ptr <io::output> setup_output (std::shared_ptr <io::output> output, int state = real_real, int flags = 0x00) {
			// Iterate through the scalar fields and append them to the variables which the output will write to file
			TRACE ("Setting up output stream...");

			if (!(params ["output.output"].as <bool> ())) return NULL;

			if (!(flags & no_variables)) {
				for (data::iterator iter = begin (); iter != end (); ++iter) {
					TRACE ("Appending " << *iter << " to output...");
					output->append (*iter, (*this) (*iter, state));
				}
				output->full=1;
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

		template <class format>
		std::shared_ptr <io::output> setup_output_from (YAML::Node output_params, const formats::data_grid &grid, int state = real_real, int flags = 0x00) {
			// Iterate through the scalar fields and append them to the variables for which the input will search
			TRACE ("Setting up output from YAML Node...");
			if (file_from (output_params) == "") return NULL;

			std::shared_ptr <io::output> output_stream;
			if (output_params ["timed"].IsDefined () && output_params ["timed"].as <bool> ()) {
				output_stream.reset (new io::timed_appender_output <datatype, formats::netcdf> (grid, file_from (output_params), duration, output_params ["timed_every"].as <datatype> ()));
			} else {
				output_stream.reset (new io::appender_output <formats::netcdf> (grid, file_from (output_params), output_params ["every"].as <int> ()));
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
			output_ptr->template append <datatype> ("t", &duration, formats::scalar);
			output_ptr->template append <datatype> ("dt", &timestep, formats::scalar);
			output_ptr->template append <const int> ("mode", &(get_mode ()), formats::scalar);
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
		
	protected:
		/*!**********************************************************************
		 * \brief Initialize a new dataset in the data object
		 ************************************************************************/
		virtual grids::variable <datatype> &_initialize (std::string i_name, int i_flags = 0x00) = 0;
	};
	
	/*!**********************************************************************
	 * \brief An implemeted form of data in 2D
	 ************************************************************************/
	template <class datatype>
	class implemented_data : public data <datatype>
	{
	protected:
		using data <datatype>::iterator;
		using data <datatype>::weights;
		using data <datatype>::variables;
		using data <datatype>::flags;
		
		std::shared_ptr <grids::grid <datatype>> grid_n; //!< The horizontal grid object
		std::shared_ptr <grids::grid <datatype>> grid_m; //!< The vertical grid object
		int n; //!< The horizontal extent of the data
		int m; //!< The vertical extent of the data

		std::vector <datatype> area;
		
	public:
		using data <datatype>::transformers;
		using data <datatype>::duration;
		
		/*!**********************************************************************
		 * \param i_axis_n The horizontal axis object
		 * \param i_axis_m The vertical axis object
		 * \param i_name The integer name of the element
		 * \param dump_file The string name of the dump file (no extension)
		 * \param dump_directory The string directory to store the dump
		 * \param dump_every The frequency to dump to file
		 ************************************************************************/
		implemented_data (grids::axis *i_axis_n, grids::axis *i_axis_m, io::parameters &i_params, int i_name = 0, std::string dump_file = "", std::string dump_directory = "./", int dump_every = 1) : 
		data <datatype> (i_params, i_name), 
		grid_n (std::shared_ptr <grids::grid <datatype>> (new typename grids::horizontal::grid <datatype> (i_axis_n))), 
		grid_m (std::shared_ptr <grids::grid <datatype>> (new typename grids::vertical::grid <datatype> (i_axis_m))), 
		n (grid_n->get_n ()), 
		m (grid_m->get_n ()) {
			// Set up output
			const formats::data_grid o_grid = formats::data_grid::two_d (n, m, 0, 0, 0, 0);
			
			std::shared_ptr <io::output> dump_stream;
			if (dump_file != "") {
				std::string file_format = dump_directory + dump_file;
				char buffer [file_format.size () * 2];
				snprintf (buffer, file_format.size () * 2, file_format.c_str (), i_name);

				dump_stream.reset (new io::replace_output <formats::netcdf> (o_grid, buffer, dump_every));
				this->setup_dump (dump_stream);
			}
			
			this->initialize ("x");
			this->initialize ("z");

			// For weighted averages, calculate area
			area.resize (n * m);
			for (int i = 1; i < n; ++i) {
				for (int j = 1; j < m; ++j) {
					area [i * m + j] = ((*(this->grid_n)) [i] - (*(this->grid_n)) [i - 1]) * ((*(this->grid_m)) [j] - (*(this->grid_m)) [j - 1]);
				}
			}

			weights = &area [0];
			
			/*
				TODO Clean dump file generation
			*/
		}
		
		virtual ~implemented_data () {}

		datatype *operator() (const std::string &name, int state = 0, int i = 0, int j = 0) {
			return data <datatype>::operator() (name, state, i * m + j);
		}
		
		/*!**********************************************************************
		 * \copydoc data::setup_profile
		 ************************************************************************/
		virtual void setup_profile (std::shared_ptr <io::output> output_ptr, int flags = 0x00) {
			typedef typename std::map <std::string, std::shared_ptr <grids::variable <datatype>>>::iterator iterator;
			for (iterator iter = data <datatype>::variables.begin (); iter != data <datatype>::variables.end (); ++iter) {
				output_ptr->template append <datatype> (iter->first, (*this) (iter->first));
				output_ptr->template append <datatype> ("rms_" + iter->first, std::shared_ptr <functors::functor> (new functors::root_mean_square_functor <datatype> ((*this) (iter->first), n, m)));
				output_ptr->template append <datatype> ("avg_" + iter->first, std::shared_ptr <functors::functor> (new functors::average_functor <datatype> ((*this) (iter->first), n, m)));
			}
			output_ptr->template append <datatype> ("t", &(duration), formats::scalar);
			output_ptr->template append <const int> ("mode", &(data <datatype>::get_mode ()), formats::scalar);

			data <datatype>::setup_profile (output_ptr, flags);
		}

		template <class format>
		std::shared_ptr <io::input> setup_from (YAML::Node input_params) {
			const formats::data_grid grid = formats::data_grid::two_d (n, m);

			return data <datatype>::template setup_from <format> (input_params, grid);
		}

		template <class format>
		std::shared_ptr <io::output> setup_output_from (YAML::Node output_params, int state = real_real, int flags = 0x00) {
			const formats::data_grid grid = formats::data_grid::two_d (n, m);

			return data <datatype>::template setup_output_from <format> (output_params, grid, state, flags);
		}

		std::shared_ptr <functors::functor> output_max (std::string variable) {
			return std::shared_ptr <functors::functor> (new functors::max_functor <double> (n, m, (*this) (variable)));
		}

		std::shared_ptr <functors::functor> output_avg (std::string variable) {
			return std::shared_ptr <functors::functor> (new functors::weighted_average_functor <double> (n, m, this->weights, (*this) (variable)));
		}

		std::shared_ptr <functors::functor> output_deriv (std::string variable) {
			return typename std::shared_ptr <functors::functor> (new functors::average_functor <double> (n, 1, typename std::shared_ptr <functors::functor> (new typename functors::slice_functor <double> (n, m, m / 2, typename std::shared_ptr <functors::functor> (new functors::deriv_functor <double> ((*this) (variable), n, m, &(*grid_m) [0]))))));
		}

		std::shared_ptr <functors::functor> output_flux (std::string variable, std::string velocity) {
			return typename std::shared_ptr <functors::functor> (new functors::average_functor <double> (n, 1, typename std::shared_ptr <functors::functor> (new typename functors::slice_functor <double> (n, m, m / 2, typename std::shared_ptr <functors::functor> (new typename functors::product_functor <double> (n, m, (*this) (velocity), (*this) (variable)))))));
		}
	
	protected:
		/*!**********************************************************************
		 * \copydoc data::_initialize
		 ************************************************************************/
		virtual grids::variable <datatype> &_initialize (std::string name, int i_flags = 0x00) {
			TRACE ("Initializing " << name << "...");
			// Size allowing for real FFT buffer
			if (i_flags & uniform_n) {
				data <datatype>::variables [name] = std::shared_ptr <grids::variable <datatype>> (new grids::variable <datatype> (*grid_m, flags, name));
				return *variables [name];
			} else if (i_flags & uniform_m) {
				data <datatype>::variables [name] = std::shared_ptr <grids::variable <datatype>> (new grids::variable <datatype> (*grid_n, flags, name));
				return *variables [name];
			}
			data <datatype>::variables [name] = std::shared_ptr <grids::variable <datatype>> (new grids::variable <datatype> (*grid_n, *grid_m, flags, name, 3, i_flags & vector2D ? 2 : (i_flags & vector3D ? 3 : 1)));
			if (name == "x") {
				for (int j = 0; j < m; ++j) {
					linalg::copy (n, &((*grid_n) [0]), (*this) (name, real_real, 0, j), 1, m);
				}
				transformers [name] = NULL;
			} else if (name == "z") {
				for (int i = 0; i < n; ++i) {
					linalg::copy (m, &((*grid_m) [0]), (*this) (name, real_real, i));
				}
				transformers [name] = NULL;
			} else {
				transformers [name] = std::shared_ptr <plans::transforms::transformer <datatype> > (new plans::transforms::implemented_transformer <datatype> (*grid_n, *grid_m, (*this) [name], plans::transforms::forward_vertical | plans::transforms::forward_horizontal | plans::transforms::inverse_vertical | plans::transforms::inverse_horizontal , &(flags), &((*this) [name].component_flags)));
			}
			TRACE ("Done.");
			return *variables [name];
		}
	};
} /* data */

#endif /* end of include guard: DATA_HPP_574F2B9E */
