/*!**********************************************************************
 * \file virtual.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef VIRTUAL_HPP_7BF7D5B0
#define VIRTUAL_HPP_7BF7D5B0

#include <map>
#include <cstring>
#include <typeinfo>
#include <array>

#include "logger/logger.hpp"

#include "format.hpp"
#include "exceptions.hpp"

namespace formats
{
	class virtual_file;

   /*!**********************************************************************
    * \brief A map of virtual files to be used like disk output, i.e. every unique string input maps to exactly one virtual_file object
    ************************************************************************/
   extern std::map <std::string, formats::virtual_file> virtual_files;
} /* formats */

namespace formats
{
	/*!**********************************************************************
	 * \brief A class designed to act as a virtual file
	 ************************************************************************/
	class virtual_file
	{
		std::map <std::string, void *> data_map; //!< A map containing void pointers to the data
		std::map <std::string, const std::type_info *> types; //!< A map containing the types of the pointers in the data_map
		std::map <std::string, size_t> sizes; //!< A map containing the sizes of the types in the data_map, for indexing
	public:
		std::map <std::string, std::array <int, 2>> dims; //!< A map containing the dimensions of the data in data_map

		/*
			TODO The dims member is a bit messy
		*/

		virtual_file () {}

		virtual ~virtual_file () {
			for (std::map <std::string, void *>::iterator iter = data_map.begin (); iter != data_map.end (); iter++) {
				free (iter->second);
			}
		}

		/*!**********************************************************************
		 * \brief We redefine the assignment operator in order to perform a deep copy
		 ************************************************************************/
		virtual_file &operator= (virtual_file &virtual_file) {
			if (&virtual_file == this) {
				return *this;
			}

			for (std::map <std::string, void *>::iterator iter = virtual_file.begin (); iter != virtual_file.end (); ++iter) {
				_add_var (iter->first, virtual_file.types [iter->first], virtual_file.sizes [iter->first], virtual_file.dims [iter->first] [0], virtual_file.dims [iter->first] [1]);
				memcpy (data_map [iter->first], iter->second, virtual_file.sizes [iter->first] * virtual_file.dims [iter->first] [0] * virtual_file.dims [iter->first] [1]);
			}

			return *this;
		}

		/*!**********************************************************************
		 * \brief Make an iterator at the start of the virtual_file keys, for iterating through the keys
		 * 
		 * I'm not certain this is a useful ability, considering that the contents can have varying types
		 ************************************************************************/
		std::map <std::string, void *>::iterator begin () {
			return data_map.begin ();
		}

		/*!**********************************************************************
		 * \brief Make an iterator at the end of the virtual_file keys
		 ************************************************************************/
		std::map <std::string, void *>::iterator end () {
			return data_map.end ();
		}

		/*!**********************************************************************
		 * \brief Reset the virtual_file object and free its contents
		 ************************************************************************/
		virtual void reset () {
			for (std::map <std::string, void *>::iterator iter = data_map.begin (); iter != data_map.end (); iter++) {
				free (iter->second);
			}
			data_map.clear ();
			types.clear ();
			sizes.clear ();
			dims.clear ();
		}

		/*!**********************************************************************
		 * \brief Get the value of the given variable for the given index
		 * 
		 * \param name The string representation of the variable
		 * \param i The integer horizontal index
		 * \param j The integer vertical index
		 * 
		 * Beware: this method is very slow, and should only be used for debugging and error handling purposes. For use of the object in actual code, use the put and get methods for optimal speed.
		 * 
		 * \return A reference to the variable value at the given index
		 ************************************************************************/
		template <class datatype>
		datatype &index (std::string name, int i = 0, int j = 0) {
			if (check_type <datatype> (name)) {
				return ((datatype *) data_map [name]) [i * dims [name] [1] + j];
			} else {
				ERROR ("Incorrect type");
				throw exceptions::bad_type ();
			}
		}

		/*!**********************************************************************
		 * \brief Add a variable to the class
		 * 
		 * \param name The string representation of the variable
		 * \param n The integer size in the horizontal dimension
		 * \param m The integer size in the vertical dimension
		 ************************************************************************/
		template <class datatype>
		void add_var (std::string name, int n = 1, int m = 1) {
			_add_var (name, &typeid (datatype), sizeof (datatype), n, m);
			datatype *temp = (datatype *) (data_map [name]);
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < m; ++j) {
					temp [i * m + j] = 0.0;
				}
			}
		}

		/*!**********************************************************************
		 * \brief Check the type of the given variable
		 * 
		 * \param name The string representation of the variable
		 * 
		 * This is mainly for error checking.
		 ************************************************************************/
		template <class datatype>
		bool check_type (std::string name) {
			if (typeid (datatype) == *types [name]) {
				return true;
			} else {
				return false;
			}
		}

		/*!**********************************************************************
		 * \brief Put data into the virtual_file object
		 * 
		 * \param name The string representation of the variable
		 * \param data The pointer to the data to copy
		 * \param n The integer horizontal size of the data array
		 * \param m The integer vertical size of the data array
		 * \param ldm The integer leading dimension size of the input array (if -1, use m)
		 * 
		 * This method currently only copies contiguous memory to contiguous memory. It could be rewritten in the future to work on subarrays. 
		 ************************************************************************/
		template <class datatype>
		void put (std::string name, datatype *data, int n = 1, int m = 1, int ldm = -1) {
			TRACE ("Putting " << name << "...");
			if (check_type <datatype> (name)) {
				ldm = ldm == -1 ? m : ldm;
				memcpy (data_map [name], data, sizeof (datatype) * n * ldm);
				TRACE ("Done.");
			} else {
				ERROR ("Incorrect type");
				throw exceptions::bad_type ();
			}
		}

		/*!**********************************************************************
		 * \brief Get the data from the virtual_file object
		 * 
		 * \param name The string representation of the variable
		 * \param data The pointer toward the final data location
		 * \param n The integer horizontal size of the data array
		 * \param m The integer vertical size of the data array
		 * \param ldm The integer leading dimension size of the input array (if -1, use m)
		 * 
		 * This method currently only copies contiguous memory to contiguous memory. It could be rewritten in the future to work on subarrays.
		 ************************************************************************/
		template <class datatype>
		void get (std::string name, datatype *data, int n = 1, int m = 1, int ldm = -1) {
			TRACE ("Getting " << name << "...");
			if (check_type <datatype> (name)) {
				ldm = ldm == -1 ? m : ldm;
				memcpy (data, data_map [name], sizeof (datatype) * n * ldm);
			} else {
				ERROR ("Incorrect type");
				throw exceptions::bad_type ();
			}
		}

	private:
		/*!**********************************************************************
		 * \brief This private member contains the implementation of add_var
		 * 
		 * \param name The string representation of the variable
		 * \param type A pointer to the std::type_info object of the given type
		 * \param size The size_t size of the given type
		 * \param n The integer horizontal size of the variable
		 * \param m The integer vertical size of the variable
		 ************************************************************************/
		void _add_var (std::string name, const std::type_info *type, size_t size, int n = 1, int m = 1) {
			TRACE ("Adding variable " << name << "...");
			if (data_map [name]) {
				free (data_map [name]);
			}
			data_map [name] = malloc (size * n * m);
			types [name] = type;
			sizes [name] = size;
			dims [name] [0] = n;
			dims [name] [1] = m;
			TRACE ("Done. " << data_map [name]);
		}
	};
	
	/*!**********************************************************************
	 * \brief A format designed to hold data in memory
	 * 
	 * This class contains only static members, so an instantiation of this class is not necessary
	 ************************************************************************/
	class virtual_format
	{
	public:
		static bool uses_files; //!< A boolean flag indicating whether this format uses files (false)
		
		virtual_format () {}
		
		~virtual_format () {}
		
		/*!**********************************************************************
		 * \return The extension associated with virtual files ("")
		 ************************************************************************/
		static std::string extension () {return "";}
		
		static bool is_open (std::string file_name) {
			if (virtual_files.find(file_name) != virtual_files.end()) return true;
			return false;
		}

		/*!**********************************************************************
		 * \brief Open the virtual file for read/write
		 ************************************************************************/
		static void open_file (const data_grid &grid, std::string file_name, int file_type) {
			// if (file_type == read_file && (!virtual_files [file_name])) {
				// ERROR ("Virtual file doesn't exist.");
				// throw 0;
			// } 
			/*
				TODO Check for virtual file existence
			*/
			virtual_files [file_name];
		}

		static void set_time (std::string file_name, double time) {}
		
		/*!**********************************************************************
		 * \brief Close the file
		 * 
		 * Nothing to do here.
		 ************************************************************************/
		static void close_file (std::string file_name, int file_type) {
		}
		
		static void add_global_attribute (std::string file_name, std::string name, std::string &attribute) {}

		/*!**********************************************************************
		 * \brief Write to file
		 * 
		 * \param grid The data_grid object containing the information on how to output
		 * \param file_name The file name to write to
		 * \param name The name of the data object to write
		 * \param data The pointer to the data to output
		 * \param record The record stamp at which to output
		 * \param flags An integer flag describing which dimensions to output
		 * 
		 * This is one of the main methods of the class, describing how to get the data into the output object
		 ************************************************************************/
		template <class datatype>
		static void write (const data_grid &grid, std::string file_name, std::string name, void *data, int record = -1, int flags = all_d) {
			int n = grid.get_n (0), m = grid.get_n (1);
			if (! (flags & one_d)) {
				n = 1;
			}
			if (! (flags & m_profile)) {
				m = 1;
			}
			virtual_files [file_name].add_var <datatype> (name, n, m);
			virtual_files [file_name].put <datatype> (name, (datatype *) data, n, m);
		}

		/*!**********************************************************************
		 * \brief Write to file
		 * 
		 * \param grid The data_grid object containing the information on how to input
		 * \param file_name The file name to write to
		 * \param name The name of the data object to write
		 * \param data The pointer to the data to output
		 * \param record The record stamp from which to input
		 * \param flags An integer flag describing which dimensions to input
		 * 
		 * This is one of the main methods of the class, describing how to get the data from the input object
		 ************************************************************************/
		template <class datatype>
		static void read (const data_grid &grid, std::string file_name, std::string name, void *data, int record = -1, int flags = all_d) {
			int n = grid.get_n (0), m = grid.get_n (1);
			if (! (flags & one_d)) {
				n = 1;
			}
			if (! (flags & m_profile)) {
				m = 1;
			}
			virtual_files [file_name].get <datatype> (name, (datatype *) data, n, m);
		}
	};
} /* formats */

#endif /* end of include guard: VIRTUAL_HPP_7BF7D5B0 */
