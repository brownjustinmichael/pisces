/*!**********************************************************************
 * \file output.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef OUTPUT_HPP_31603039
#define OUTPUT_HPP_31603039

#include <fstream>

#include "logger/logger.hpp"

#include "formats/format.hpp"
#include "formats/exceptions.hpp"
#include "functors/functor.hpp"

#include "versions/version.hpp"

/*!*******************************************************************
 * \namespace io
 * 
 * \brief A namespace containing all the input and output classes of 
 * the code
 *********************************************************************/
namespace io
{
	/*!*******************************************************************
	 * \brief An abstract output stream base class that generates output files
	 * 
	 * There's a lot of overlap between this and input. Perhaps we should design a superclass for the two or merge them into one class.
	 *********************************************************************/
	class output
	{
	protected:
		std::string file_name; //!< The string file name
		int file_format; //!< The file format from io_flags: e.g. replace_file
		const formats::data_grid grid; //!< The data_grid object associated with the output
		
		typedef void func_t (const formats::data_grid &grid, std::string file_name, std::string name, void *data, int record, int dims); //!< A prototype for write functions
		std::vector <func_t *> write_functions; //!< A vector of pointers to the write functions
		std::vector <std::string> names; //!< A vector of the string representations of the variables
		std::vector <void *> data_ptrs; //!< A vector of pointers to the arrays of data
		std::vector <std::shared_ptr <functors::functor>> functor_ptrs; //!< A vector of pointers to the functors
		std::vector <int> dims; //!< A vector of the integer number of dimensions for each output

	public:
		/*!*******************************************************************
		 * \param i_grid The data_grid object representing the structure of the data
		 * \param i_file_name The string representation of the output file; do not include the extension; it will be added later
		 * \param i_file_format The integer io_flag associated with the desired output type (e.g. replace_file)
		 *********************************************************************/
		output (formats::data_grid i_grid, std::string i_file_name = "out", int i_file_format = formats::replace_file) : file_name (i_file_name), file_format (i_file_format), grid (i_grid) {
			DEBUG ("Names has size " << names.size ());
		}

		virtual ~output () {}
		
		/*!**********************************************************************
		 * \return The version of the class
		 ************************************************************************/
		static versions::version& version () {
			static versions::version version ("1.1.0.0");
			return version;
		}
		
		/*
			TODO The get_function formality is rather obnoxious; perhaps there's a better way
		*/
		
		/*!**********************************************************************
		 * \brief Given a float pointer, return the appropriate write_function
		 * 
		 * \return A pointer to the appropriate write function for float arrays
		 ************************************************************************/
		virtual func_t *get_function (const float *) = 0;

		/*!**********************************************************************
		 * \brief Given a double pointer, return the appropriate write_function
		 * 
		 * \return A pointer to the appropriate write function for double arrays
		 ************************************************************************/
		virtual func_t *get_function (const double *) = 0;
		
		/*!**********************************************************************
		 * \brief Given an int pointer, return the appropriate write_function
		 * 
		 * \return A pointer to the appropriate write function for int arrays
		 ************************************************************************/
		virtual func_t *get_function (const int *) = 0;

		/*!**********************************************************************
		 * \brief Append a datatype functor to the list to be output
		 * 
		 * \param name The string representation of the quantity
		 * \param functor_ptr A pointer to the functor that will calculate the resulting quantity
		 * \param flags A set of binary flags indicating which dimensions to output
		 * 
		 * WARNING: THIS WILL HAPPILY ACCEPT A FUNCTOR OF NON-DATATYPE TEMPLATE
		 ************************************************************************/
		template <class datatype>
		void append (std::string name, std::shared_ptr <functors::functor> functor_ptr, int flags = formats::all_d) {
			TRACE ("Appending " << name << " to output...");
			for (int i = 0; i < (int) names.size (); ++i) {
				if (names [i] == name) {
					WARN ("Reuse of name " << name);
					functor_ptrs [i] = functor_ptr;
					return;
				}
			}
			functor_ptrs.push_back (functor_ptr);
			append <datatype> (name, (datatype *) functor_ptr->calculate (), flags);
			TRACE ("Functor appended.");
		}

		/*!*******************************************************************
		 * \brief Append a datatype array to the list to be output
		 * 
		 * \param name The string representation of the variable
		 * \param data_ptr A datatype pointer to the data
		 * \param flags A set of flags indicating what dimensions to output, defaults to all available dimensions
		 * 
		 * Since this is shared between input and output, we might give them a shared superclass in the future.
		 *********************************************************************/
		template <class datatype>
		void append (std::string name, datatype *data_ptr, int flags = formats::all_d) {
			TRACE ("Appending " << name << " to output..." << *data_ptr);
			for (int i = 0; i < (int) names.size (); ++i) {
				if (names [i] == name) {
					WARN ("Reuse of name " << name);
					data_ptrs [i] = (void *) data_ptr;
					return;
				}
			}
			names.push_back (name);
			dims.push_back (flags);
			data_ptrs.push_back ((void *) data_ptr);
			write_functions.push_back (get_function (data_ptr));
			TRACE ("Appended.");
		}
		
		/*!*******************************************************************
		 * \brief A function to output the data to file
		 * 
		 * This function should be overwritten by subclasses.
		 *********************************************************************/
		virtual void to_file (int record = -1) = 0;
	};
	
	/*!**********************************************************************
	 * \brief An output class that takes a format object as a template argument
	 * 
	 * Note that the format object should be modeled after the ones in the formats folder. They require the static template methods extension, open_file, write, and close_file. The choice of this somewhat circumlocuting format is to effectively allow the formats to be of the same abstract type but have different implementations of the same template functions. (Since abstract template functions don't make much sense in C++)
	 * 
	 * \copydoc output
	 ************************************************************************/
	template <class format>
	class formatted_output : public output
	{
	public:
		/*!**********************************************************************
		 * \copydoc output::output
		 ************************************************************************/
		formatted_output (formats::data_grid i_grid, std::string i_file_name = "out", int i_file_format = formats::replace_file) : output (i_grid, i_file_name + format::extension (), i_file_format) {}		

		virtual ~formatted_output () {}
		
		/*!**********************************************************************
		 * \copydoc output::get_function
		 * 
		 * This finds the appropriate write function in the format template argument and returns it.
		 ************************************************************************/
		virtual func_t *get_function (const float *ptr) {
			return &format::template write <float>;
		}
		
		/*!**********************************************************************
		 * \copydoc output::get_function
		 * 
		 * This finds the appropriate write function in the format template argument and returns it.
		 ************************************************************************/
		virtual func_t *get_function (const double *ptr) {
			return &format::template write <double>;
		}
		
		/*!**********************************************************************
		 * \copydoc output::get_function
		 * 
		 * This finds the appropriate write function in the format template argument and returns it.
		 ************************************************************************/
		virtual func_t *get_function (const int *ptr) {
			return &format::template write <int>;
		}
		
		/*!**********************************************************************
		 * \brief Check it the file opens; otherwise, raise an exception
		 * 
		 * This is in the output class to avoid having to check if files open correctly in the format classes.
		 ************************************************************************/
		virtual void check_file (std::string file_name) {
			std::ifstream filestr;
			filestr.open (file_name);
			if (!(filestr.is_open ())) {
				std::ofstream filestr;
				filestr.open (file_name);
				if (!(filestr.is_open ())) {
					throw formats::exceptions::file_exception (file_name);
				}
				filestr.close ();
			}
			filestr.close ();
		}
		
		/*!**********************************************************************
		 * \copybrief output::to_file
		 ************************************************************************/
		virtual void to_file (int record = -1) {
			TRACE ("Sending to file...");
	
			INFO ("Outputting to file " << file_name << "...");
			
			if (format::uses_files) check_file (file_name.c_str ());
			format::open_file (grid, file_name.c_str (), output::file_format);
	
			// Calculate the inner values of any relevant functors
			for (int i = 0; i < (int) functor_ptrs.size (); ++i) {
				functor_ptrs [i]->calculate ();
			}
			
			// Output the array data
			for (int i = 0; i < (int) write_functions.size (); ++i) {
				write_functions [i] (grid, file_name, names [i], data_ptrs [i], record, dims [i]);
			}
	
			format::close_file (file_name.c_str (), output::file_format);
		}
	};

	/*!**********************************************************************
	 * \brief An output class that increments file names each output
	 * 
	 * \copydoc formatted_output
	 ************************************************************************/
	template <class format>
	class incremental : public formatted_output <format>
	{
	private:
		std::string file_format; //!< The string format to create the file names
		int output_every; //!< The integer frequency of outputs
		int count; //!< The current integer count of total outputs

	public:
		/*!**********************************************************************
		 * \param i_grid The data_grid object representing the structure of the data
		 * \param i_file_format A string file format, using standard string formatting for the increments (e.g. "output_%02d")
		 * \param i_output_every The integer frequency of outputs
		 ************************************************************************/
		incremental (formats::data_grid i_grid, std::string i_file_format, int i_output_every = 1) : formatted_output <format> (i_grid, "", formats::replace_file), file_format (i_file_format + format::extension ()), output_every (i_output_every > 0 ? i_output_every : 1),
		count (0) {}

		virtual ~incremental () {}

		/*!**********************************************************************
		 * \copybrief formatted_output::to_file
		 * 
		 * Also increments the file count and checks whether to output at all
		 ************************************************************************/
		void to_file (int record = -1) {
			TRACE ("Sending to file...");
			if (count % output_every == 0) {
				char buffer [file_format.size () * 2];
				snprintf (buffer, file_format.size () * 2, file_format.c_str (), count / output_every);
				formatted_output <format>::file_name = buffer;
				formatted_output <format>::to_file (record);
			}
			++count;
		}
	};

	/*!**********************************************************************
	 * \brief An output class that appends each output to the same file
	 * 
	 * \copydoc formatted_output
	 ************************************************************************/
	template <class format>
	class appender_output : public formatted_output <format>
	{
	private:
		int output_every; //!< The integer frequency of outputs
		int count; //!< The current integer count of total outputs

	public:
		/*!**********************************************************************
		 * \param i_grid The data_grid object representing the structure of the data
		 * \param i_file_name A string file name
		 * \param i_output_every The integer frequency of outputs
		 ************************************************************************/
		appender_output (formats::data_grid i_grid, std::string i_file_name, int i_output_every = 1) : formatted_output <format> (i_grid, i_file_name, formats::append_file), output_every (i_output_every > 0 ? i_output_every : 1), count (0) {
			DEBUG ("CONSTRUCTING");
		}

		virtual ~appender_output () {}

		/*!**********************************************************************
		 * \copybrief formatted_output::to_file
		 * 
		 * Also increments the file count and checks whether to output at all
		 ************************************************************************/
		void to_file (int record = -1) {
			TRACE ("Sending to file...");
			if (count % output_every == 0) {
				formatted_output <format>::to_file (record ? record : count);
			}
			++count;
		}
	};
	
	/*!**********************************************************************
	 * \brief An output class that appends each output to the same file
	 * 
	 * \copydoc formatted_output
	 ************************************************************************/
	template <class datatype, class format>
	class timed_appender_output : public formatted_output <format>
	{
	private:
		datatype &duration;
		datatype output_every; //!< The integer frequency of outputs
		datatype previous;

	public:
		/*!**********************************************************************
		 * \param i_grid The data_grid object representing the structure of the data
		 * \param i_file_name A string file name
		 * \param i_output_every The integer frequency of outputs
		 ************************************************************************/
		timed_appender_output (formats::data_grid i_grid, std::string i_file_name, datatype &i_duration, datatype i_output_every = 1.0) : formatted_output <format> (i_grid, i_file_name, formats::append_file), duration (i_duration), output_every (i_output_every > 0.0 ? i_output_every : 1.0), previous (-output_every) {}

		virtual ~timed_appender_output () {}

		/*!**********************************************************************
		 * \copybrief formatted_output::to_file
		 * 
		 * Also increments the file count and checks whether to output at all
		 ************************************************************************/
		void to_file (int record = -1) {
			TRACE ("Sending to file...");
			while (duration - previous >= output_every) {
				formatted_output <format>::to_file (record);
				previous += output_every;
			}
		}
	};
	
	/*!**********************************************************************
	 * \brief An output class that replaces the previous output
	 ************************************************************************/
	template <class format>
	class replace_output : public formatted_output <format>
	{
	private:
		int output_every; //!< The integer frequency of outputs
		int count; //!< The current integer count of total outputs

	public:
		/*!**********************************************************************
		 * \param i_grid The data_grid object representing the structure of the data
		 * \param i_file_name The string representation of the output file; do not include the extension; it will be added later
		 * \param i_output_every The integer frequency of outputs
		 ************************************************************************/
		replace_output (formats::data_grid i_grid, std::string i_file_name, int i_output_every = 1) : formatted_output <format> (i_grid, i_file_name, formats::replace_file), output_every (i_output_every > 0 ? i_output_every : 1), count (0) {}

		virtual ~replace_output () {}

		/*!**********************************************************************
		 * \copybrief formatted_output::to_file
		 * 
		 * Also increments the file count and checks whether to output at all
		 ************************************************************************/
		void to_file (int record = -1) {
			TRACE ("Sending to file...");
			if (count % output_every == 0) {
				formatted_output <format>::to_file ();
			}
			++count;
		}
	};
} /* io */

#endif /* end of include guard: OUTPUT_HPP_31603039 */
