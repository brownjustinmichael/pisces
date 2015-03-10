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

#include "versions/version.hpp"
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
		formats::data_grid grid; //!< The data_grid object that keeps track of the numerical structure of the data
		std::vector <std::string> names; //!< A vector of the string representations of the variables
		typedef void func_t (const formats::data_grid &grid, std::string file_name, std::string name, void *data, int record, int dims); //!< A typedef of the read_function prototype
		std::vector <func_t *> read_functions; //!< A vector containing all the functions to execute upon read
		std::vector <int> dims; //!< A vector containing the relevant dimensions of the read evaluations
		std::vector <void *> data_ptrs; //!< A vector of pointers to the arrays of data
	
	public:
		/*!**********************************************************************
		 * /copydoc output::output
		 ************************************************************************/
		input (formats::data_grid i_grid, std::string i_file_name = "in") : file_name (i_file_name), grid (i_grid) {}
		
		virtual ~input () {}
		
		/*!**********************************************************************
		 * \brief Get the version of the class
		 ************************************************************************/
		static versions::version& version () {
			static versions::version version ("1.1.0.0");
			return version;
		}
		
		/*!**********************************************************************
		 * \brief Given a float pointer, return the appropriate read_function
		 * 
		 * \return A pointer to the appropriate read function for float arrays
		 ************************************************************************/
		virtual func_t *get_function (const float *) = 0;
		
		/*!**********************************************************************
		 * \brief Given a double pointer, return the appropriate read_function
		 * 
		 * \return A pointer to the appropriate read function for double arrays
		 ************************************************************************/
		virtual func_t *get_function (const double *) = 0;
		
		/*!**********************************************************************
		 * \brief Given a float pointer, return the appropriate read_function
		 * 
		 * \return A pointer to the appropriate read function for int arrays
		 ************************************************************************/
		virtual func_t *get_function (const int *) = 0;

		/*!*******************************************************************
		 * \brief Append a datatype array to the list to be output
		 * 
		 * \param name The string representation of the variable to append
		 * \param data_ptr A datatype pointer to the data array
		 *********************************************************************/
		template <class datatype>
		void append (std::string name, datatype *data_ptr, int flags = formats::all_d) {
			TRACE ("Appending " << name << " to input...");
			for (int i = 0; i < (int) names.size (); ++i) {
				if (names [i] == name) {
					WARN ("Reuse of name " << name);
					data_ptrs [i] = (void *) data_ptr;
					return;
				}
			}
			names.push_back (name);
			data_ptrs.push_back ((void *) data_ptr);
			dims.push_back (flags);
			read_functions.push_back (get_function (data_ptr));
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
	 * Note that the format object should be modeled after the ones in this file. They require the static methods extension, write, open_file, and close_file.
	 ************************************************************************/
	template <class format>
	class formatted_input : public input
	{
	public:
		/*!**********************************************************************
		 * \copydoc input::input
		 ************************************************************************/
		formatted_input (formats::data_grid i_grid, std::string i_file_name = "in") :
		input (i_grid, i_file_name + format::extension ()) {}
	
		virtual ~formatted_input () {}
		
		/*!**********************************************************************
		 * \copydoc input::get_function
		 * 
		 * This finds the appropriate read function in the format template argument and returns it.
		 ************************************************************************/
		virtual func_t *get_function (const float *ptr) {
			return &format::template read <float>;
		}
		
		/*!**********************************************************************
		 * \copydoc input::get_function
		 * 
		 * This finds the appropriate read function in the format template argument and returns it.
		 ************************************************************************/
		virtual func_t *get_function (const double *ptr) {
			return &format::template read <double>;
		}
		
		/*!**********************************************************************
		 * \copydoc input::get_function
		 * 
		 * This finds the appropriate read function in the format template argument and returns it.
		 ************************************************************************/
		virtual func_t *get_function (const int *ptr) {
			return &format::template read <int>;
		}
	
		/*!**********************************************************************
		 * \copybrief input::from_file
		 ************************************************************************/
		virtual void from_file (int record = -1) {
			INFO ("Inputting from file " << file_name << "...");
		
			format::open_file (grid, file_name.c_str (), formats::read_file);
			
			for (int i = 0; i < (int) read_functions.size (); ++i) {
				read_functions [i] (grid, file_name, names [i], data_ptrs [i], record, dims [i]);
			}
		
			format::close_file (file_name.c_str (), formats::read_file);
		}
	};
} /* io */

#endif /* end of include guard: INPUT_HPP_EFD01D95 */
