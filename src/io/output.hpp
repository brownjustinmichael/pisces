/*!**********************************************************************
 * \file output.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef OUTPUT_HPP_31603039
#define OUTPUT_HPP_31603039

#include "logger/logger.hpp"

#include "formats/format.hpp"
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
		std::vector <std::string> names; //!< A vector of the string representations of the variables
		std::vector <std::string> scalar_names; //!< A vector of the string representations of the scalar variables
		std::vector <std::string> functor_names; //!< A vector of the string representations of the functor variables
		std::vector <std::string> scalar_functor_names; //!< A vector of the string representations of the functor variables
		std::vector <const std::type_info*> types; //!< A vector of the types of the variables
		std::vector <const std::type_info*> scalar_types; //!< A vector of the types of the scalar variables
		std::vector <const std::type_info*> functor_types; //!< A vector of the types of the functor variables
		std::vector <const std::type_info*> scalar_functor_types; //!< A vector of the types of the functor variables
		std::vector <void *> data_ptrs; //!< A vector of pointers to the arrays of data
		std::vector <void *> scalar_ptrs; //!< A vector of pointers to the scalar data
		std::vector <void *> functor_ptrs; //!< A vector of pointers to the functors
		std::vector <void *> scalar_functor_ptrs; //!< A vector of pointers to the functors

	public:
		/*!*******************************************************************
		 * \param i_file_name The string representation of the output file; do not include the extension; it will be added later
		 * 
		 * This seems an extremely tedious way to do this. It might be superior to accept arrays or brace initialized lists
		 *********************************************************************/
		output (data_grid i_grid, std::string i_file_name = "out", int i_file_format = replace_file) : file_name (i_file_name), file_format (i_file_format), grid (i_grid) {}

		virtual ~output () {}

		/*!*******************************************************************
		 * \brief Append a datatype array to the list to be output
		 * 
		 * \param name The string representation of the variable
		 * \param data_ptr A datatype pointer to the data
		 * 
		 * Since this is shared between input and output, we might give them a shared superclass in the future.
		 *********************************************************************/
		template <class datatype>
		void append (std::string name, datatype *data_ptr) {
			TRACE ("Appending " << name << " to output..." << *data_ptr);
			for (int i = 0; i < (int) names.size (); ++i) {
				if (names [i] == name) {
					WARN ("Reuse of name " << name);
					types [i] = &typeid (datatype);
					data_ptrs [i] = (void *) data_ptr;
					return;
				}
			}
			names.push_back (name);
			types.push_back (&typeid (datatype));
			data_ptrs.push_back ((void *) data_ptr);
			TRACE ("Appended.");
		}

		/*!**********************************************************************
		 * \brief Append a datatype scalar to the list to be output
		 * 
		 * \param name The string representation of the variable
		 * \param data_ptr A datatype pointer to the data
		 * 
		 * Since this is shared between input and output, we might give them a shared superclass in the future.
		 ************************************************************************/
		template <class datatype>
		void append_scalar (std::string name, datatype *data_ptr) {
			TRACE ("Appending " << name << " to output...");
			for (int i = 0; i < (int) scalar_names.size (); ++i) {
				if (scalar_names [i] == name) {
					WARN ("Reuse of name " << name);
					scalar_types [i] = &typeid (datatype);
					scalar_ptrs [i] = (void *) data_ptr;
					TRACE ("Scalar updated.");
					return;
				}
			}
			scalar_types.push_back (&typeid (datatype));
			scalar_names.push_back (name);
			scalar_ptrs.push_back ((void *) data_ptr);
			TRACE ("Scalar appended.");
		}

		/*!**********************************************************************
		 * \brief Append a datatype functor to the list to be output
		 * 
		 * \param name The string representation of the quantity
		 * \param functor_ptr A pointer to the functor that will calculate the resulting quantity
		 ************************************************************************/
		template <class datatype>
		void append_functor (std::string name, functors::functor <datatype> *functor_ptr) {
			TRACE ("Appending " << name << " to output...");
			for (int i = 0; i < (int) functor_names.size (); ++i) {
				if (functor_names [i] == name) {
					WARN ("Reuse of name " << name);
					functor_types [i] = &typeid (datatype);
					functor_ptrs [i] = (void *) functor_ptr;
					TRACE ("Scalar updated.");
					return;
				}
			}
			functor_types.push_back (&typeid (datatype));
			functor_names.push_back (name);
			functor_ptrs.push_back ((void *) functor_ptr);
			TRACE ("Functor appended.");
		}

		/*!**********************************************************************
		 * \brief Append a datatype scalar_functor to the list to be output
		 * 
		 * \param name The string representation of the quantity
		 * \param scalar_functor_ptr A pointer to the scalar_functor that will calculate the resulting quantity
		 ************************************************************************/
		template <class datatype>
		void append_scalar_functor (std::string name, functors::functor <datatype> *scalar_functor_ptr) {
			TRACE ("Appending " << name << " to output...");
			for (int i = 0; i < (int) scalar_functor_names.size (); ++i) {
				if (scalar_functor_names [i] == name) {
					WARN ("Reuse of name " << name);
					scalar_functor_types [i] = &typeid (datatype);
					scalar_functor_ptrs [i] = (void *) scalar_functor_ptr;
					TRACE ("Scalar updated.");
					return;
				}
			}
			scalar_functor_types.push_back (&typeid (datatype));
			scalar_functor_names.push_back (name);
			scalar_functor_ptrs.push_back ((void *) scalar_functor_ptr);
			TRACE ("Functor appended.");
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

		/*!**********************************************************************
		 * \copybrief output::to_file
		 ************************************************************************/
		virtual void to_file (int record = -1) {
			TRACE ("Sending to file...");
	
			INFO ("Outputting to file " << file_name << "...");
	
			format::open_file (grid, file_name.c_str (), output::file_format);
	
			// Output the scalars
			for (int i = 0; i < (int) scalar_names.size (); ++i) {
				if (scalar_types [i] == &typeid (double)) {
					format::template write_scalar <double> (file_name, scalar_names [i], (double *) scalar_ptrs [i], record);
				} else if (scalar_types [i] == &typeid (float)) {
					format::template write_scalar <float> (file_name, scalar_names [i], (float *) scalar_ptrs [i], record);
				} else if (scalar_types [i] == &typeid (int)) {
					format::template write_scalar <int> (file_name, scalar_names [i], (int *) scalar_ptrs [i], record);
				} else {
					throw 0;
				}
			}
	
			// Output the scalar_functors
			for (int i = 0; i < (int) scalar_functor_names.size (); ++i) {
				if (scalar_functor_types [i] == &typeid (double)) {
					format::template write_scalar <double> (file_name, scalar_functor_names [i], ((functors::functor <double> *) scalar_functor_ptrs [i])->calculate (), record);
				} else if (scalar_functor_types [i] == &typeid (float)) {
					format::template write_scalar <float> (file_name, scalar_functor_names [i], ((functors::functor <float> *) scalar_functor_ptrs [i])->calculate (), record);
				} else if (scalar_functor_types [i] == &typeid (int)) {
					format::template write_scalar <int> (file_name, scalar_functor_names [i], ((functors::functor <int> *) scalar_functor_ptrs [i])->calculate (), record);
				} else {
					throw 0;
				}
			}
	
			// Output the array data
			for (int i = 0; i < (int) names.size (); ++i) {
				if (types [i] == &typeid (double)) {
					format::template write <double> (grid, file_name, names [i], (double *) data_ptrs [i], record);
				} else if (types [i] == &typeid (float)) {
					format::template write <float> (grid, file_name, names [i], (float *) data_ptrs [i], record);
				} else if (types [i] == &typeid (int)) {
					format::template write <int> (grid, file_name, names [i], (int *) data_ptrs [i], record);
				} else {
					throw 0;
				}
			}
	
			// Output the results from the functors
			for (int i = 0; i < (int) functor_names.size (); ++i) {
				if (functor_types [i] == &typeid (double)) {
					format::template write <double> (grid, file_name, functor_names [i], ((functors::functor <double> *) functor_ptrs [i])->calculate (), record);
				} else if (functor_types [i] == &typeid (float)) {
					format::template write <float> (grid, file_name, functor_names [i], ((functors::functor <float> *) functor_ptrs [i])->calculate (), record);
				} else if (functor_types [i] == &typeid (int)) {
					format::template write <int> (grid, file_name, functor_names [i], ((functors::functor <int> *) functor_ptrs [i])->calculate (), record);
				} else {
					throw 0;
				}
			}
			if (output::file_format != append_file) {
				format::close_file (file_name.c_str (), output::file_format);
			}
	
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
} /* io */

#endif /* end of include guard: OUTPUT_HPP_31603039 */
