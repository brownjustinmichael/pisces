/*!***********************************************************************
 * \file io.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef IO_HPP_C1E9B6EF
#define IO_HPP_C1E9B6EF

#include <vector>
#include <array>
#include <map>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <memory>
#include <typeinfo>
#include "exceptions.hpp"
#include "../config.hpp"
#include <yaml-cpp/yaml.h>
#include <netcdf>
#include "utils.hpp"

namespace io
{	
	/*!**********************************************************************
	 * \brief A set of io flags to be used with the input/output classes specifying the file type
	 ************************************************************************/
	enum io_flags {
		read_file = 0,
		replace_file = 1,
		append_file = 2
	};

	class virtual_dump;
	
	/*!**********************************************************************
	 * \brief A map of virtual dumps to be used like disk output, i.e. every unique string input maps to exactly one virtual_dump object
	 ************************************************************************/
	extern std::map <std::string, virtual_dump> virtual_dumps;
	
	/*!**********************************************************************
	 * \brief Returns the netCDF object associated with a given type
	 * 
	 * \param type A pointer to a type_info object associated with the object type in question
	 * 
	 * This function throws an io::bad_type exception if the type is not known
	 * 
	 * \return The netCDF::NcType object associated with the input type
	 ************************************************************************/
	inline netCDF::NcType netcdf_type (const std::type_info* type) {
		if (type == &typeid (double)) {
			return netCDF::ncDouble;
		} else if (type == &typeid (int)) {
			return netCDF::ncInt;
		} else if (type == &typeid (float)) {
			return netCDF::ncFloat;
		} else {
			FATAL ("Unrecognized NetCDF type");
			throw exceptions::bad_type ();
		}
	}
	
	/*!**********************************************************************
	 * \brief A class designed to act as a virtual dump file
	 ************************************************************************/
	class virtual_dump
	{
		std::map <std::string, void *> data_map; //!< A map containing void pointers to the data
		std::map <std::string, const std::type_info *> types; //!< A map containing the types of the pointers in the data_map
		std::map <std::string, size_t> sizes; //!< A map containing the sizes of the types in the data_map, for indexing
	public:
		std::map <std::string, std::array <int, 2>> dims; //!< A map containing the dimensions of the data in data_map
		
		/*
			TODO The dims member is a bit messy
		*/
		
		virtual_dump () {}
		
		virtual ~virtual_dump () {
			for (std::map <std::string, void *>::iterator iter = data_map.begin (); iter != data_map.end (); iter++) {
				free (iter->second);
			}
		}
		
		/*!**********************************************************************
		 * \brief We redefine the assignment operator in order to perform a deep copy
		 ************************************************************************/
		virtual_dump &operator= (virtual_dump &dump) {
			if (&dump == this) {
				return *this;
			}
		
			for (std::map <std::string, void *>::iterator iter = dump.begin (); iter != dump.end (); ++iter) {
				_add_var (iter->first, dump.types [iter->first], dump.sizes [iter->first], dump.dims [iter->first] [0], dump.dims [iter->first] [1]);
				memcpy (data_map [iter->first], iter->second, dump.sizes [iter->first] * dump.dims [iter->first] [0] * dump.dims [iter->first] [1]);
			}
		
			return *this;
		}
		
		/*!**********************************************************************
		 * \brief Make an iterator at the start of the virtual_dump keys, for iterating through the keys
		 * 
		 * I'm not certain this is a useful ability, considering that the contents can have varying types
		 ************************************************************************/
		std::map <std::string, void *>::iterator begin () {
			return data_map.begin ();
		}
		
		/*!**********************************************************************
		 * \brief Make an iterator at the end of the virtual_dump keys
		 ************************************************************************/
		std::map <std::string, void *>::iterator end () {
			return data_map.end ();
		}
		
		/*!**********************************************************************
		 * \brief Reset the virtual_dump object and free its contents
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
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < m; ++j) {
					((datatype *) (data_map [name])) [i * m + j] = 0.0;
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
		 * \brief Put data into the virtual_dump object
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
		 * \brief Get the data from the virtual_dump object
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
				TRACE ("Done.");
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
	 * \brief An extension of the YAML::Node class for convenience
	 * 
	 * This implementation uses the convenient shorthand of the '.' character delimiting nested YAML trees. This means that this character should be avoided in parameter names in your YAML file.
	 ************************************************************************/
	class parameters : public YAML::Node
	{
	public:
		using YAML::Node::operator=; //!< Use the assignment operator from the super class
		
		/*!**********************************************************************
		 * \param file_name The parameter file from which the parameters should be loaded
		 ************************************************************************/
		parameters (std::string file_name = "") {
			if (file_name != "") {
				YAML::Node::operator= (YAML::LoadFile (file_name));
			}
		}
		
		virtual ~parameters () {}
		
		/*!**********************************************************************
		 * \brief Overload the index operator for the YAML::Node
		 * 
		 * \param key The string representation of the parameter
		 * 
		 * Note that this returns a copy of a YAML::Node object, which can be used to check and set parameters, but these objects do not treat parameter delimiters like this class.
		 * 
		 * \return A copy of the YAML::Node object at the given parameter.
		 ************************************************************************/
		YAML::Node operator[] (std::string key);
		
		/*!**********************************************************************
		 * \brief Get the parameter associated with the given key
		 * 
		 * \param key The string representation of the parameter
		 * 
		 * This method is included for convenience. It avoids some of the messy YAML syntax and takes advantage of the overwritted index operator. This method throws a key_does_not_exist exception if it cannot find the key.
		 * 
		 * \return Returns the value of the given parameter
		 ************************************************************************/
		template <typename datatype>
		const datatype get (std::string key) {
			if (operator[] (key).IsDefined ()) {
				return operator[] (key).as <datatype> ();
			} else {
				throw exceptions::key_does_not_exist (key);
			}
		}
	};
	
	/*!**********************************************************************
	 * \brief Abstract class for the format_functor object
	 * 
	 * The functor class is designed to take an instruction and apply it to some data for visualization and optimization purposes. For example, a functor could take a two dimensional grid of data and produce a one dimensional average or profile. This class serves as a wrapper for the calculate function, which returns a pointer to the processed data.
	 ************************************************************************/
	template <class datatype>
	class format_functor
	{
	public:
		virtual ~format_functor () {}
		
		/*!**********************************************************************
		 * \brief The instruction to process on the data
		 * 
		 * The class serves as a wrapper for this function
		 * 
		 * \return A pointer to the processed data, for output
		 ************************************************************************/
		virtual datatype *calculate () = 0;
	};
	
	/*!**********************************************************************
	 * \brief Averages a two dimensional block of data
	 ************************************************************************/
	template <class datatype>
	class weighted_average_functor : public format_functor <datatype>
	{
	private:
		datatype *weight, *data; //!< A datatype pointer to the input data
		format_functor <datatype> *functor;
		int n; //!< The integer horizontal extent of the data
		int m; //!< The integer vertical extent of the data
		datatype inner_data; //!< A vector of processed data to output
	
	public:
		/*!**********************************************************************
		 * \param i_data The datatype pointer to the data to average
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 ************************************************************************/
		weighted_average_functor (int i_n, int i_m, datatype *i_weight, datatype *i_data) : weight (i_weight), data (i_data), n (i_n), m (i_m) {
		}
	
		weighted_average_functor (int i_n, int i_m, datatype *i_weight, format_functor <datatype> *i_functor) : weight (i_weight), data (i_functor->calculate ()), functor (i_functor), n (i_n), m (i_m) {
		}
	
		/*!**********************************************************************
		 * \brief Average the data and return a pointer to the first element
		 * 
		 * \return The first element of the averaged 1D array
		 ************************************************************************/
		datatype *calculate () {
			if (functor) {
				functor->calculate ();
			}
			inner_data = (datatype) 0;
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < n; ++i) {
					inner_data += weight [i * m + j] * data [i * m + j];
				}
			}
			return &inner_data;
		}
	};
	
	/*!**********************************************************************
	 * \brief Averages a two dimensional block of data in the horizontal direction
	 ************************************************************************/
	template <class datatype>
	class average_functor : public format_functor <datatype>
	{
	private:
		datatype *data; //!< A datatype pointer to the input data
		std::shared_ptr <format_functor <datatype>> functor;
		int n; //!< The integer horizontal extent of the data
		int m; //!< The integer vertical extent of the data
		std::vector <datatype> inner_data; //!< A vector of processed data to output
		
	public:
		/*!**********************************************************************
		 * \param i_data The datatype pointer to the data to average
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 ************************************************************************/
		average_functor (datatype *i_data, int i_n, int i_m) : data (i_data), n (i_n), m (i_m) {
			inner_data.resize (m);
		}
		
		average_functor (format_functor <datatype> *i_functor, int i_n, int i_m) : data (i_functor->calculate ()), functor (std::shared_ptr <format_functor <datatype>> (i_functor)), n (i_n), m (i_m) {
			inner_data.resize (m);
		}
		
		/*!**********************************************************************
		 * \brief Average the data and return a pointer to the first element
		 * 
		 * \return The first element of the averaged 1D array
		 ************************************************************************/
		datatype *calculate () {
			if (functor) {
				functor->calculate ();
			}
			for (int j = 0; j < m; ++j) {
				inner_data [j] = (datatype) 0;
				for (int i = 0; i < n; ++i) {
					inner_data [j] += data [i * m + j];
				}
				inner_data [j] /= (datatype) n;
			}
			return &inner_data [0];
		}
	};
	
	/*!**********************************************************************
	 * \brief Finds the root-mean-square of data in the horizontal direction
	 ************************************************************************/
	template <class datatype>
	class root_mean_square_functor : public format_functor <datatype>
	{
	private:
		datatype *data; //!< A datatype pointer to the input data
		int n; //!< The integer horizontal extent of the data
		int m; //!< The integer vertical extent of the data
		std::vector <datatype> inner_data; //!< A datatype vector of processed data to output
		
	public:
		/*!**********************************************************************
		 * \param i_data The datatype pointer to the data to root-mean-square
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 ************************************************************************/
		root_mean_square_functor (datatype *i_data, int i_n, int i_m) : data (i_data), n (i_n), m (i_m) {
			inner_data.resize (m);
		}
		
		/*!**********************************************************************
		 * \brief Take the root-mean-square and return a pointer to the first element
		 * 
		 * \return The first element of the resulting array
		 ************************************************************************/
		datatype *calculate () {
			for (int j = 0; j < m; ++j) {
				inner_data [j] = (datatype) 0;
				for (int i = 0; i < n; ++i) {
					inner_data [j] += data [i * m + j] * data [i * m + j];
				}
				inner_data [j] = sqrt(inner_data [j] / (datatype) n);
			}
			return &inner_data [0];
		}
	};
	
	/*!**********************************************************************
	 * \brief Averages a two dimensional block of data in the horizontal direction
	 ************************************************************************/
	template <class datatype>
	class div_functor : public format_functor <datatype>
	{
	private:
		datatype *pos_x, *pos_z, *data_x, *data_z; //!< A datatype pointer to the input data
		int n; //!< The integer horizontal extent of the data
		int m; //!< The integer vertical extent of the data
		std::vector <datatype> inner_data; //!< A vector of processed data to output
	
	public:
		/*!**********************************************************************
		 * \param i_data The datatype pointer to the data to average
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 ************************************************************************/
		div_functor (datatype *i_pos_x, datatype *i_pos_z, datatype *i_data_x, datatype *i_data_z, int i_n, int i_m) : pos_x (i_pos_x), pos_z (i_pos_z), data_x (i_data_x), data_z (i_data_z), n (i_n), m (i_m) {
			inner_data.resize (n * m);
		}
	
		/*!**********************************************************************
		 * \brief Average the data and return a pointer to the first element
		 * 
		 * \return The first element of the averaged 1D array
		 ************************************************************************/
		datatype *calculate () {
			for (int i = 0; i < n; ++i) {
				inner_data [i * m] = (data_z [i * m + 1] - data_z [i * m]) / (pos_z [i * m + 1] - pos_z [i * m]);
				for (int j = 1; j < m - 1; ++j) {
					inner_data [i * m + j] = (data_z [i * m + j + 1] - data_z [i * m + j - 1]) / (pos_z [i * m + j + 1] - pos_z [i * m + j - 1]);
				}
				inner_data [i * m + m - 1] = (data_z [i * m + m - 1] - data_z [i * m + m - 2]) / (pos_z [i * m + m - 1] - pos_z [i * m + m - 2]);
			}
			for (int j = 0; j < m; ++j) {
				inner_data [j] += (data_x [m + j] - data_x [(n - 1) * m + j]) / (pos_x [m + j] - pos_x [(n - 1) * m + j]);
				for (int i = 1; i < n - 1; ++i) {
					inner_data [i * m + j] += (data_x [(i + 1) * m + j] - data_x [(i - 1) * m + j]) / (pos_x [(i + 1) * m + j] - pos_x [(i - 1) * m + j]);
				}
				inner_data [(n - 1) * m + j] += (data_x [j] - data_x [(n - 2) * m + j]) / (pos_x [j] - pos_x [(n - 2) * m + j]);
			}
			return &inner_data [0];
		}
	};
	
	/*!**********************************************************************
	 * \brief Averages a two dimensional block of data in the horizontal direction
	 ************************************************************************/
	template <class datatype>
	class product_functor : public format_functor <datatype>
	{
	private:
		datatype *data_1, *data_2; //!< A datatype pointer to the input data
		int n; //!< The integer horizontal extent of the data
		int m; //!< The integer vertical extent of the data
		std::vector <datatype> inner_data; //!< A vector of processed data to output

	public:
		/*!**********************************************************************
		 * \param i_data The datatype pointer to the data to average
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 ************************************************************************/
		product_functor (int i_n, int i_m, datatype *i_data_1, datatype *i_data_2) : data_1 (i_data_1), data_2 (i_data_2), n (i_n), m (i_m) {
			inner_data.resize (n * m);
		}

		/*!**********************************************************************
		 * \brief Average the data and return a pointer to the first element
		 * 
		 * \return The first element of the averaged 1D array
		 ************************************************************************/
		datatype *calculate () {
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < m; ++j) {
					inner_data [i * m + j] = data_1 [i * m + j] * data_2 [i * m + j];
				}
			}
			return &inner_data [0];
		}
	};
	
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
		 * \param i_n The integer number of points in the first dimension of the data
		 * \param i_m The integer number of points in the second dimension of the data
		 * \param i_l The integer number of points in the third dimension of the data
		 * \param i_n_max The integer total extent of the first dimension (including fill space; if 0, use i_n)
		 * \param i_m_max The integer total extent of the second dimension (including fill space; if 0, use i_m)
		 * \param i_l_max The integer total extent of the third dimension (including fill space; if 0, use i_l)
		 * \param i_n_offset The integer offset of the starting index in the first dimension
		 * \param i_m_offset The integer offset of the starting index in the second dimension
		 * \param i_l_offset The integer offset of the starting index in the third dimension
		 * 
		 * This seems an extremely tedious way to do this. It might be superior to accept arrays or brace initialized lists
		 *********************************************************************/
		output (std::string i_file_name = "out", int i_file_format = replace_file, int i_n = 1, int i_m = 1, int i_l = 1, int i_n_max = 0, int i_m_max = 0, int i_l_max = 0, int i_n_offset = 0, int i_m_offset = 0, int i_l_offset = 0) : file_name (i_file_name), file_format (i_file_format), n (i_n), m (i_m), l (i_l), n_max (i_n_max ? i_n_max : n), m_max (i_m_max ? i_m_max : m), l_max (i_l_max ? i_l_max : l), n_offset (i_n_offset), m_offset (i_m_offset), l_offset (i_l_offset) {}
		
		/*
			TODO This is a mess... the abstract output shouldn't need to know about the dimensions
		*/
		
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
		void append_functor (std::string name, format_functor <datatype> *functor_ptr) {
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
		void append_scalar_functor (std::string name, format_functor <datatype> *scalar_functor_ptr) {
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
		formatted_output (std::string i_file_name = "out", int i_file_format = replace_file, int i_n = 1, int i_m = 1, int i_l = 1, int i_n_max = 0, int i_m_max = 0, int i_l_max = 0, int i_n_offset = 0, int i_m_offset = 0, int i_l_offset = 0) : output (i_file_name + format::extension (), i_file_format, i_n, i_m, i_l, i_n_max, i_m_max, i_l_max, i_n_offset, i_m_offset, i_l_offset) {}		
		
		virtual ~formatted_output () {
			format::close_file (file_name.c_str (), output::file_format);
		}
		
		/*!**********************************************************************
		 * \copybrief output::to_file
		 ************************************************************************/
		virtual void to_file (int record = -1) {
			TRACE ("Sending to file...");
			
			INFO ("Outputting to file " << file_name << "...");
			
			format::open_file (file_name.c_str (), output::file_format, n_max, m_max, l_max);
			
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
					format::template write_scalar <double> (file_name, scalar_functor_names [i], ((format_functor <double> *) scalar_functor_ptrs [i])->calculate (), record);
				} else if (scalar_functor_types [i] == &typeid (float)) {
					format::template write_scalar <float> (file_name, scalar_functor_names [i], ((format_functor <float> *) scalar_functor_ptrs [i])->calculate (), record);
				} else if (scalar_functor_types [i] == &typeid (int)) {
					format::template write_scalar <int> (file_name, scalar_functor_names [i], ((format_functor <int> *) scalar_functor_ptrs [i])->calculate (), record);
				} else {
					throw 0;
				}
			}
			
			// Output the array data
			for (int i = 0; i < (int) names.size (); ++i) {
				if (types [i] == &typeid (double)) {
					format::template write <double> (file_name, names [i], (double *) data_ptrs [i], n, m, l, n_offset, m_offset, l_offset, record);
				} else if (types [i] == &typeid (float)) {
					format::template write <float> (file_name, names [i], (float *) data_ptrs [i], n, m, l, n_offset, m_offset, l_offset, record);
				} else if (types [i] == &typeid (int)) {
					format::template write <int> (file_name, names [i], (int *) data_ptrs [i], n, m, l, n_offset, m_offset, l_offset, record);
				} else {
					throw 0;
				}
			}
			
			// Output the results from the functors
			for (int i = 0; i < (int) functor_names.size (); ++i) {
				if (functor_types [i] == &typeid (double)) {
					format::template write <double> (file_name, functor_names [i], ((format_functor <double> *) functor_ptrs [i])->calculate (), n, m, l, n_offset, m_offset, l_offset, record);
				} else if (functor_types [i] == &typeid (float)) {
					format::template write <float> (file_name, functor_names [i], ((format_functor <float> *) functor_ptrs [i])->calculate (), n, m, l, n_offset, m_offset, l_offset, record);
				} else if (functor_types [i] == &typeid (int)) {
					format::template write <int> (file_name, functor_names [i], ((format_functor <int> *) functor_ptrs [i])->calculate (), n, m, l, n_offset, m_offset, l_offset, record);
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
		incremental (std::string i_file_format, int i_output_every = 1, int i_n = 1, int i_m = 1, int i_l = 1, int i_n_max = 0, int i_m_max = 0, int i_l_max = 0, int i_n_offset = 0, int i_m_offset = 0, int i_l_offset = 0) :
		formatted_output <format> ("", replace_file, i_n, i_m, i_l, i_n_max, i_m_max, i_l_max, i_n_offset, i_m_offset, i_l_offset),
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
		appender_output (std::string i_file_name, int i_output_every = 1, int i_n = 1, int i_m = 1, int i_l = 1, int i_n_max = 0, int i_m_max = 0, int i_l_max = 0, int i_n_offset = 0, int i_m_offset = 0, int i_l_offset = 0) : formatted_output <format> (i_file_name, append_file, i_n, i_m, i_l, i_n_max, i_m_max, i_l_max, i_n_offset, i_m_offset, i_l_offset), output_every (i_output_every > 0 ? i_output_every : 1), count (0) {}
	
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

#endif /* end of include guard: IO_HPP_C1E9B6EF */
