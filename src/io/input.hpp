/*!**********************************************************************
 * \file input.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef INPUT_HPP_EFD01D95
#define INPUT_HPP_EFD01D95

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
		int n; //!< The integer number of points in the first dimension of the data
		int m; //!< The integer number of points in the second dimension of the data
		int l; //!< The integer number of points in the third dimension of the data
		int n_max; //!< The integer total extent of the first dimension of the data
		int m_max; //!< The integer total extent of the second dimension of the data
		int l_max; //!< The integer total extent of the third dimension of the data
		int n_offset; //!< The integer offset of the starting index in the first dimension
		int m_offset; //!< The integer offset of the starting index in the second dimension
		int l_offset; //!< The integer offset of the starting index in the third dimension
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
		input (std::string i_file_name = "in", int i_n = 1, int i_m = 1, int i_l = 1, int i_n_max = 0, int i_m_max = 0, int i_l_max = 0, int i_n_offset = 0, int i_m_offset = 0, int i_l_offset = 0) : file_name (i_file_name), n (i_n), m (i_m), l (i_l), n_max (i_n_max), m_max (i_m_max), l_max (i_l_max), n_offset (i_n_offset), m_offset (i_m_offset), l_offset (i_l_offset) {}
	
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
		formatted_input (std::string i_file_name = "in", int i_n = 1, int i_m = 1, int i_l = 1, int i_n_max = 0, int i_m_max = 0, int i_l_max = 0, int i_n_offset = 0, int i_m_offset = 0, int i_l_offset = 0) :
		input (i_file_name + format::extension (), i_n, i_m, i_l, i_n_max, i_m_max, i_l_max, i_n_offset, i_m_offset, i_l_offset) {}
	
		virtual ~formatted_input () {}
	
		/*!**********************************************************************
		 * \copybrief input::from_file
		 ************************************************************************/
		virtual void from_file (int record = -1) {
			INFO ("Inputting from file " << file_name << "...");
		
			format::open_file (file_name.c_str (), read_file, n_max, m_max, l_max);
		
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
					format::template read <double> (file_name, names [i], (double *) data_ptrs [i], n, m, l, n_offset, m_offset, l_offset, record);
				} else if (types [i] == &typeid (float)) {
					format::template read <float> (file_name, names [i], (float *) data_ptrs [i], n, m, l, n_offset, m_offset, l_offset, record);
				} else if (types [i] == &typeid (int)) {
					format::template read <int> (file_name, names [i], (int *) data_ptrs [i], n, m, l, n_offset, m_offset, l_offset, record);
				} else {
					throw 0;
				}
			}
		
			format::close_file (file_name.c_str (), read_file);
		}
	};
} /* io */

#endif /* end of include guard: INPUT_HPP_EFD01D95 */
