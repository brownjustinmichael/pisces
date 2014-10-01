/*!**********************************************************************
 * \file input.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef INPUT_HPP_EFD01D95
#define INPUT_HPP_EFD01D95

#include "logger/logger.hpp"

#include "formats/format.hpp"
#include "functors/functor.hpp"

namespace io
{
	/*!**********************************************************************
	 * \brief An abstract base class for file input
	 * 
	 * There's a lot of overlap between this and output. Perhaps we should design a superclass for the two or merge them into one class.
	 ************************************************************************/
	class input
	{
	protected:
		std::string file_name; //!< The string representation of the file
		data_grid grid;
		std::vector <std::string> names; //!< A vector of the string representations of the variables
		std::vector <std::string> scalar_names; //!< A vector of the string representations of the scalar variables
		std::vector <const std::type_info*> types; //!< A vector of the types of the variables
		std::vector <const std::type_info*> scalar_types; //!< A vector of the types of the scalar variables
		std::vector <void *> data_ptrs; //!< A vector of pointers to the arrays of data
		std::vector <void *> scalar_ptrs; //!< A vector of pointers to the scalar data
	
	public:
		/*!**********************************************************************
		 * /copydoc output::output
		 ************************************************************************/
		input (data_grid i_grid, std::string i_file_name = "in") : file_name (i_file_name), grid (i_grid) {}
	
		virtual ~input () {}
	
		/*!*******************************************************************
		 * \brief Append a datatype array to the list to be output
		 * 
		 * \param name The string representation of the variable to append
		 * \param data_ptr A datatype pointer to the data array
		 *********************************************************************/
		template <class datatype>
		void append (std::string name, datatype *data_ptr) {
			TRACE ("Appending " << name << " to input...");
			for (int i = 0; i < (int) names.size (); ++i) {
				if (names [i] == name) {
					WARN ("Reuse of name " << name);
					types [i] = &typeid (datatype);
					data_ptrs [i] = (void *) data_ptr;
					return;
				}
			}
			types.push_back (&typeid (datatype));
			names.push_back (name);
			data_ptrs.push_back ((void *) data_ptr);
		}
	
		/*!*******************************************************************
		 * \brief Append a datatype scalar to the list to be output
		 * 
		 * \param name The string representation of the variable to append
		 * \param data_ptr A datatype pointer to the data value
		 *********************************************************************/
		template <class datatype>
		void append_scalar (std::string name, datatype *data_ptr) {
			TRACE ("Appending " << name << " to input...");
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
		}
	
		/*!*******************************************************************
		 * \brief A virtual function to output the data to file
		 * 
		 * This function should be overwritten by subclasses.
		 *********************************************************************/
		virtual void from_file (int record = -1) = 0;
	};
	
	/*!**********************************************************************
	 * \brief An input class that takes a format object as a template argument
	 * 
	 * Note that the format object should be modeled after the ones in this file. They require the static methods extension, write, write_scalar, and write_functor.
	 ************************************************************************/
	template <class format>
	class formatted_input : public input
	{
	public:
		/*!**********************************************************************
		 * \copydoc input::input
		 ************************************************************************/
		formatted_input (data_grid i_grid, std::string i_file_name = "in") :
		input (i_grid, i_file_name + format::extension ()) {}
	
		virtual ~formatted_input () {}
	
		/*!**********************************************************************
		 * \copybrief input::from_file
		 ************************************************************************/
		virtual void from_file (int record = -1) {
			INFO ("Inputting from file " << file_name << "...");
		
			format::open_file (grid, file_name.c_str (), read_file);
		
			// Input the scalars from file
			for (int i = 0; i < (int) scalar_names.size (); ++i) {
				if (scalar_types [i] == &typeid (double)) {
					format::template read_scalar <double> (file_name, scalar_names [i], (double *) scalar_ptrs [i], record);
				} else if (scalar_types [i] == &typeid (float)) {
					format::template read_scalar <float> (file_name, scalar_names [i], (float *) scalar_ptrs [i], record);
				} else if (scalar_types [i] == &typeid (int)) {
					format::template read_scalar <int> (file_name, scalar_names [i], (int *) scalar_ptrs [i], record);
				} else {
					throw 0;
				}
			}

			// Input the arrays from file
			for (int i = 0; i < (int) names.size (); ++i) {
				if (types [i] == &typeid (double)) {
					format::template read <double> (grid, file_name, names [i], (double *) data_ptrs [i], record);
				} else if (types [i] == &typeid (float)) {
					format::template read <float> (grid, file_name, names [i], (float *) data_ptrs [i], record);
				} else if (types [i] == &typeid (int)) {
					format::template read <int> (grid, file_name, names [i], (int *) data_ptrs [i], record);
				} else {
					throw 0;
				}
			}
		
			format::close_file (file_name.c_str (), read_file);
		}
	};
} /* io */

#endif /* end of include guard: INPUT_HPP_EFD01D95 */
