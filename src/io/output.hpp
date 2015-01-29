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
		int file_format;
		const data_grid grid;
		typedef void func_t (const data_grid &grid, std::string file_name, std::string name, void *data, int record, int dims);
		std::vector <func_t *> write_functions;
		std::vector <std::string> names; //!< A vector of the string representations of the variables
		std::vector <void *> data_ptrs; //!< A vector of pointers to the arrays of data
		std::vector <std::shared_ptr <functors::functor>> functor_ptrs; //!< A vector of pointers to the functors
		std::vector <int> dims;

	public:
		/*!*******************************************************************
		 * \param i_file_name The string representation of the output file; do not include the extension; it will be added later
		 * 
		 * This seems an extremely tedious way to do this. It might be superior to accept arrays or brace initialized lists
		 *********************************************************************/
		output (data_grid i_grid, std::string i_file_name = "out", int i_file_format = replace_file) : file_name (i_file_name), file_format (i_file_format), grid (i_grid) {}

		virtual ~output () {}
		
		virtual func_t *get_function (const float *) = 0;
		virtual func_t *get_function (const double *) = 0;
		virtual func_t *get_function (const int *) = 0;

		/*!**********************************************************************
		 * \brief Append a datatype functor to the list to be output
		 * 
		 * \param name The string representation of the quantity
		 * \param functor_ptr A pointer to the functor that will calculate the resulting quantity
		 * 
		 * WARNING: THIS WILL HAPPILY ACCEPT A FUNCTOR OF NON-DATATYPE TEMPLATE
		 ************************************************************************/
		template <class datatype>
		void append (std::string name, std::shared_ptr <functors::functor> functor_ptr, int flags = all_d) {
			TRACE ("Appending " << name << " to output...");
			for (int i = 0; i < (int) names.size (); ++i) {
				if (names [i] == name) {
					WARN ("Reuse of name " << name);
					functor_ptrs [i] = functor_ptr;
					return;
				}
			}
			DEBUG ("Pointer " << &*functor_ptr);
			functor_ptrs.push_back (functor_ptr);
			append <datatype> (name, (datatype *) functor_ptr->calculate (), flags);
			TRACE ("Functor appended.");
		}

		/*!*******************************************************************
		 * \brief Append a datatype array to the list to be output
		 * 
		 * \param name The string representation of the variable
		 * \param data_ptr A datatype pointer to the data
		 * 
		 * Since this is shared between input and output, we might give them a shared superclass in the future.
		 *********************************************************************/
		template <class datatype>
		void append (std::string name, datatype *data_ptr, int flags = all_d) {
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
	 * Note that the format object should be modeled after the ones in this file. They require the static methods extension, write, write_scalar, and write_functor.
	 ************************************************************************/
	template <class format>
	class formatted_output : public output
	{
	public:
		/*!**********************************************************************
		 * \copydoc output::output
		 ************************************************************************/
		formatted_output (data_grid i_grid, std::string i_file_name = "out", int i_file_format = replace_file) : output (i_grid, i_file_name + format::extension (), i_file_format) {}		

		virtual ~formatted_output () {
			format::close_file (file_name.c_str (), output::file_format);
		}
		
		virtual func_t *get_function (const float *ptr) {
			return &format::template write <float>;
		}
		virtual func_t *get_function (const double *ptr) {
			return &format::template write <double>;
		}
		virtual func_t *get_function (const int *ptr) {
			return &format::template write <int>;
		}
		
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
	
			// Output the scalar_functors
			for (int i = 0; i < (int) functor_ptrs.size (); ++i) {
				functor_ptrs [i]->calculate ();
			}
			
			// Output the array data
			for (int i = 0; i < (int) write_functions.size (); ++i) {
				write_functions [i] (grid, file_name, names [i], data_ptrs [i], record, dims [i]);
			}
	
			format::close_file (file_name.c_str (), output::file_format);
	
			/*
				TODO This behavior means that in a crash, all output data are lost, appender files should be opened and closed like all others
			*/
		}
	};

	/*!**********************************************************************
	 * \brief An output class that increments file names each output
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
		 * \param i_file_format A string file format, using standard string formatting for the increments (e.g. "output_%02d")
		 * \param i_output_every The integer frequency of outputs
		 * \copydoc formatted_output::formatted_output
		 ************************************************************************/
		incremental (data_grid i_grid, std::string i_file_format, int i_output_every = 1) :
		formatted_output <format> (i_grid, "", replace_file),
		file_format (i_file_format + format::extension ()),
		output_every (i_output_every > 0 ? i_output_every : 1),
		count (0) {}

		virtual ~incremental () {}

		/*!**********************************************************************
		 * \copybrief formatted_output::to_file
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
	 ************************************************************************/
	template <class format>
	class appender_output : public formatted_output <format>
	{
	private:
		int output_every; //!< The integer frequency of outputs
		int count; //!< The current integer count of total outputs

	public:
		/*!**********************************************************************
		 * \param i_output_every The integer frequency of outputs
		 * \copydoc formatted_output::formatted_output
		 ************************************************************************/
		appender_output (data_grid i_grid, std::string i_file_name, int i_output_every = 1) : formatted_output <format> (i_grid, i_file_name, append_file), output_every (i_output_every > 0 ? i_output_every : 1), count (0) {}

		virtual ~appender_output () {}

		/*!**********************************************************************
		 * \copybrief formatted_output::to_file
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
		 * \param i_output_every The integer frequency of outputs
		 * \copydoc formatted_output::formatted_output
		 ************************************************************************/
		replace_output (data_grid i_grid, std::string i_file_name, int i_output_every = 1) : formatted_output <format> (i_grid, i_file_name, replace_file), output_every (i_output_every > 0 ? i_output_every : 1), count (0) {}

		virtual ~replace_output () {}

		/*!**********************************************************************
		 * \copybrief formatted_output::to_file
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
