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
	enum io_flags {
		read_file = 0,
		replace_file = 1,
		append_file = 2
	};

	class virtual_dump;
	
	extern std::map <std::string, virtual_dump> virtual_dumps;
	
	inline netCDF::NcType netcdf_type (const std::type_info* type) {
		if (type == &typeid (double)) {
			return netCDF::ncDouble;
		} else if (type == &typeid (int)) {
			return netCDF::ncInt;
		} else if (type == &typeid (float)) {
			return netCDF::ncFloat;
		} else {
			FATAL ("Unrecognized MPI type");
			throw 0;
		}
	}
	
	class virtual_dump
	{
	public:
		virtual_dump () {}
		
		virtual ~virtual_dump () {
			for (std::map <std::string, void *>::iterator iter = data_map.begin (); iter != data_map.end (); iter++) {
				free (iter->second);
			}
		}
		
		virtual void reset () {
			for (std::map <std::string, void *>::iterator iter = data_map.begin (); iter != data_map.end (); iter++) {
				free (iter->second);
			}
			data_map.clear ();
			types.clear ();
			sizes.clear ();
			dims.clear ();
		}
		
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
		
		std::map <std::string, void *>::iterator begin () {
			return data_map.begin ();
		}
		
		std::map <std::string, void *>::iterator end () {
			return data_map.end ();
		}
		
		template <class datatype>
		void add_var (std::string name, int n = 1, int m = 1) {
			_add_var (name, &typeid (datatype), sizeof (datatype), n, m);
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < m; ++j) {
					((datatype *) (data_map [name])) [i * m + j] = 0.0;
				}
			}
		}
		
		template <class datatype>
		bool check_type (std::string name) {
			if (typeid (datatype) == *types [name]) {
				return true;
			} else {
				return false;
			}
		}
		
		template <class datatype>
		void put (std::string name, void *data, int n = 1, int m = 1, int ldm = -1) {
			TRACE ("Putting " << name << "...");
			if (check_type <datatype> (name)) {
				ldm = ldm == -1 ? m : ldm;
				memcpy (data_map [name], data, sizeof (datatype) * n * ldm);
				TRACE ("Done.");
			} else {
				ERROR ("Incorrect type");
				throw 0;
			}
		}
		
		template <class datatype>
		void get (std::string name, void *data, int n = 1, int m = 1, int ldm = -1) {
			TRACE ("Getting " << name << "...");
			if (check_type <datatype> (name)) {
				ldm = ldm == -1 ? m : ldm;
				memcpy (data, data_map [name], sizeof (datatype) * n * ldm);
				TRACE ("Done.");
			} else {
				ERROR ("Incorrect type");
				throw 0;
			}
		}
		
		template <class datatype>
		datatype &index (std::string name, int i = 0, int j = 0) {
			if (check_type <datatype> (name)) {
				return ((datatype *) data_map [name]) [i * dims [name] [1] + j];
			} else {
				ERROR ("Incorrect type");
				throw 0;
			}
		}
					
		std::map <std::string, std::array <int, 2>> dims;
	private:
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
		
		std::map <std::string, void *> data_map;
		std::map <std::string, const std::type_info *> types;
		std::map <std::string, size_t> sizes;
	};
	
	class parameters : public YAML::Node
	{
	public:
		parameters (std::string filename = "") {
			if (filename != "") {
				YAML::Node::operator= (YAML::LoadFile (filename));
			}
		}
		
		virtual ~parameters () {}
		
		YAML::Node operator [] (std::string key);
		
		template <typename datatype>
		const datatype get (std::string key) {
			if (operator[] (key).IsDefined ()) {
				return operator [] (key).as <datatype> ();
			} else {
				throw exceptions::key_does_not_exist (key);
			}
		}
		
		using YAML::Node::operator=;
	};
	
	
	template <class datatype>
	class format_functor
	{
	public:
		virtual ~format_functor () {}
	
		virtual datatype *calculate () = 0;
	};
	
	template <class datatype>
	class average_functor : public format_functor <datatype>
	{
	public:
		average_functor (datatype *i_data, int i_n, int i_m) :
		data (i_data),
		n (i_n),
		m (i_m) {
			inner_data.resize (m);
		}
		
		datatype *calculate () {
			for (int j = 0; j < m; ++j) {
				inner_data [j] = (datatype) 0;
				for (int i = 0; i < n; ++i) {
					inner_data [j] += data [i * m + j];
				}
				inner_data [j] /= (datatype) n;
			}
			return &inner_data [0];
		}
	
	private:
		datatype *data;
		int n, m;
		std::vector <datatype> inner_data;
	};
	
	template <class datatype>
	class root_mean_square_functor : public format_functor <datatype>
	{
	public:
		root_mean_square_functor (datatype *i_data, int i_n, int i_m) :
		data (i_data),
		n (i_n),
		m (i_m) {
			inner_data.resize (m);
		}
		
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
	
	private:
		datatype *data;
		int n, m;
		std::vector <datatype> inner_data;
	};
	
	/*!*******************************************************************
	 * \brief An abstract output stream base class that generates output files
	 *********************************************************************/
	class output
	{
	public:
		/*!*******************************************************************
		 * \param i_header_ptr A pointer to the header object
		 * \param i_n The integer number of points in the data
		 *********************************************************************/
		output (std::string i_file_name = "out", int i_n = 1, int i_m = 1, int i_l = 1, int i_n_max = 0, int i_m_max = 0, int i_l_max = 0, int i_n_offset = 0, int i_m_offset = 0, int i_l_offset = 0) :
		file_name (i_file_name),
		n (i_n), m (i_m), l (i_l), n_max (i_n_max ? i_n_max : n), m_max (i_m_max ? i_m_max : m), l_max (i_l_max ? i_l_max : l), n_offset (i_n_offset), m_offset (i_m_offset), l_offset (i_l_offset) {}
		
		virtual ~output () {}
		
		/*!*******************************************************************
		 * \brief Append a datatype array to the list to be output
		 * 
		 * \param data_ptr A datatype pointer to the data to be the new column
		 *********************************************************************/
		template <class datatype>
		void append (std::string name, datatype *data_ptr) {
			TRACE ("Appending " << name << " to output...");
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
		
		/*!*******************************************************************
		 * \brief A virtual function to output the data to file
		 * 
		 * This function should be overwritten by subclasses, though it may 
		 * contain a call to this function, which will output with default 
		 * datatype representation in C++.
		 *********************************************************************/
		virtual void to_file () = 0;
		
	protected:
		std::string file_name;
		int n, m, l;
		int n_max, m_max, l_max;
		int n_offset, m_offset, l_offset;
		std::vector <std::string> names;
		std::vector <std::string> scalar_names;
		std::vector <std::string> functor_names;
		std::vector <const std::type_info*> types;
		std::vector <const std::type_info*> scalar_types;
		std::vector <const std::type_info*> functor_types;
		std::vector <void *> data_ptrs; //!< A vector of integer pointers to the arrays of data
		std::vector <void *> scalar_ptrs; //!< A vector of integer pointers to the arrays of data
		std::vector <void *> functor_ptrs;
	};
	
	template <class format>
	class formatted_output : public output
	{
	public:
		formatted_output (std::string i_file_name = "out", int i_n = 1, int i_m = 1, int i_l = 1, int i_n_max = 0, int i_m_max = 0, int i_l_max = 0, int i_n_offset = 0, int i_m_offset = 0, int i_l_offset = 0) :
		output (i_file_name + format::extension (), i_n, i_m, i_l, i_n_max, i_m_max, i_l_max, i_n_offset, i_m_offset, i_l_offset) {}		
		
		virtual ~formatted_output () {}
		
		virtual void to_file () {
			TRACE ("Sending to file...");
			
			INFO ("Outputting to file " << file_name << "...");
			
			format::open_file (file_name.c_str (), replace_file, n_max, m_max, l_max);
			
			for (int i = 0; i < (int) scalar_names.size (); ++i) {
				if (scalar_types [i] == &typeid (double)) {
					format::template write_scalar <double> (file_name, scalar_names [i], (double *) scalar_ptrs [i]);
				} else if (scalar_types [i] == &typeid (float)) {
					format::template write_scalar <float> (file_name, scalar_names [i], (float *) scalar_ptrs [i]);
				} else if (scalar_types [i] == &typeid (int)) {
					format::template write_scalar <int> (file_name, scalar_names [i], (int *) scalar_ptrs [i]);
				} else {
					throw 0;
				}
			}

			for (int i = 0; i < (int) names.size (); ++i) {
				if (types [i] == &typeid (double)) {
					format::template write <double> (file_name, names [i], (double *) data_ptrs [i], n, m, l, n_offset, m_offset, l_offset);
				} else if (types [i] == &typeid (float)) {
					format::template write <float> (file_name, names [i], (float *) data_ptrs [i], n, m, l, n_offset, m_offset, l_offset);
				} else if (types [i] == &typeid (int)) {
					format::template write <int> (file_name, names [i], (int *) data_ptrs [i], n, m, l, n_offset, m_offset, l_offset);
				} else {
					throw 0;
				}
			}
			
			for (int i = 0; i < (int) functor_names.size (); ++i) {
				if (functor_types [i] == &typeid (double)) {
					format::template write <double> (file_name, functor_names [i], ((format_functor <double> *) functor_ptrs [i])->calculate (), n, m, l, n_offset, m_offset, l_offset);
				} else if (functor_types [i] == &typeid (float)) {
					format::template write <float> (file_name, functor_names [i], ((format_functor <float> *) functor_ptrs [i])->calculate (), n, m, l, n_offset, m_offset, l_offset);
				} else if (functor_types [i] == &typeid (int)) {
					format::template write <int> (file_name, functor_names [i], ((format_functor <int> *) functor_ptrs [i])->calculate (), n, m, l, n_offset, m_offset, l_offset);
				} else {
					throw 0;
				}
			}
			format::close_file (file_name.c_str ());
		}
	};
	
	template <class format>
	class incremental : public formatted_output <format>
	{
	public:
		incremental (std::string i_file_format, int i_output_every = 1, int i_n = 1, int i_m = 1, int i_l = 1, int i_n_max = 0, int i_m_max = 0, int i_l_max = 0, int i_n_offset = 0, int i_m_offset = 0, int i_l_offset = 0) :
		formatted_output <format> ("", i_n, i_m, i_l, i_n_max, i_m_max, i_l_max, i_n_offset, i_m_offset, i_l_offset),
		file_format (i_file_format + format::extension ()),
		output_every (i_output_every > 0 ? i_output_every : 1),
		count (0) {}
		
		virtual ~incremental () {}
	
		void to_file () {
			TRACE ("Sending to file...");
			if (count % output_every == 0) {
				char buffer [file_format.size () * 2];
				snprintf (buffer, file_format.size () * 2, file_format.c_str (), count / output_every);
				formatted_output <format>::file_name = buffer;
				formatted_output <format>::to_file ();
			}
			++count;
		}
	
	private:
		std::string file_format;
		int output_every;
		int count;
	};
	
	class input
	{
	public:
		input (std::string i_file_name = "in", int i_n = 1, int i_m = 1, int i_l = 1, int i_n_max = 1, int i_m_max = 1, int i_l_max = 1, int i_n_offset = 0, int i_m_offset = 0, int i_l_offset = 0) :
		file_name (i_file_name),
		n (i_n), m (i_m), l (i_l), n_max (i_n_max), m_max (i_m_max), l_max (i_l_max), n_offset (i_n_offset), m_offset (i_m_offset), l_offset (i_l_offset) {}
		
		virtual ~input () {
			// printf ("Destroying input\n");
		}
		
		/*!*******************************************************************
		 * \brief Append a datatype array to the list to be output
		 * 
		 * \param data_ptr A datatype pointer to the data to be the new column
		 *********************************************************************/
		template <class datatype>
		void append (std::string name, datatype *data_ptr) {
			// TRACE ("Appending " << name << " to input...");
			types.push_back (&typeid (datatype));
			names.push_back (name);
			data_ptrs.push_back ((void *) data_ptr);
		}
		
		template <class datatype>
		void append_scalar (std::string name, datatype *data_ptr) {
			// TRACE ("Appending " << name << " to output...");
			scalar_types.push_back (&typeid (datatype));
			scalar_names.push_back (name);
			scalar_ptrs.push_back ((void *) data_ptr);
		}
		
		/*!*******************************************************************
		 * \brief A virtual function to output the data to file
		 * 
		 * This function should be overwritten by subclasses, though it may 
		 * contain a call to this function, which will output with default 
		 * datatype representation in C++.
		 *********************************************************************/
		virtual void from_file () = 0;
		
	protected:
		std::string file_name;
		int n, m, l;
		int n_max, m_max, l_max;
		int n_offset, m_offset, l_offset;
		std::vector <std::string> names;
		std::vector <const std::type_info*> types;
		std::vector <void *> data_ptrs; //!< A vector of integer pointers to the arrays of data
		std::vector <std::string> scalar_names;
		std::vector <const std::type_info*> scalar_types;
		std::vector <void *> scalar_ptrs; //!< A vector of integer pointers to the arrays of data
	};
	
	template <class format>
	class formatted_input : public input
	{
	public:
		formatted_input (std::string i_file_name = "in", int i_n = 1, int i_m = 1, int i_l = 1, int i_n_max = 1, int i_m_max = 1, int i_l_max = 1, int i_n_offset = 0, int i_m_offset = 0, int i_l_offset = 0) :
		input (i_file_name + format::extension (), i_n, i_m, i_l, i_n_max, i_m_max, i_l_max, i_n_offset, i_m_offset, i_l_offset) {}
		
		virtual ~formatted_input () {}
		
		virtual void from_file () {
			INFO ("Inputting from file " << file_name << "...");

			format::open_file (file_name.c_str (), read_file, n_max, m_max, l_max);
			
			for (int i = 0; i < (int) names.size (); ++i) {
				if (types [i] == &typeid (double)) {
					format::template read <double> (file_name, names [i], (double *) data_ptrs [i], n, m, l, n_offset, m_offset, l_offset);
				} else if (types [i] == &typeid (float)) {
					format::template read <float> (file_name, names [i], (float *) data_ptrs [i], n, m, l, n_offset, m_offset, l_offset);
				} else if (types [i] == &typeid (int)) {
					format::template read <int> (file_name, names [i], (int *) data_ptrs [i], n, m, l, n_offset, m_offset, l_offset);
				} else {
					throw 0;
				}
			}
			
			for (int i = 0; i < (int) scalar_names.size (); ++i) {
				DEBUG ("THING " << scalar_names [i]);
				if (scalar_types [i] == &typeid (double)) {
					format::template read_scalar <double> (file_name, scalar_names [i], (double *) scalar_ptrs [i]);
				} else if (scalar_types [i] == &typeid (float)) {
					format::template read_scalar <float> (file_name, scalar_names [i], (float *) scalar_ptrs [i]);
				} else if (scalar_types [i] == &typeid (int)) {
					format::template read_scalar <int> (file_name, scalar_names [i], (int *) scalar_ptrs [i]);
					DEBUG ("After " << scalar_names [i] << " " << *((int *) scalar_ptrs [i]));
				} else {
					throw 0;
				}
				DEBUG ("Done");
			}
			
			format::close_file (file_name.c_str ());
		}
	};
	
	namespace one_d
	{
		// /*!*******************************************************************
		//  * \brief A simple implementation of the output class
		//  * 
		//  * This class is a simple implementation of the output class.
		//  *********************************************************************/
		// class ascii : public format <ascii>
		// {
		// public:
		// 	/*!*******************************************************************
		// 	 * \param i_file_name The string name of file for output
		// 	 * \param i_n The integer number of points in the data
		// 	 * \param i_output_every An integer number of steps between outputs
		// 	 *********************************************************************/
		// 	ascii (int i_n, std::string i_comment = "#") : 
		// 	n (i_n),
		// 	comment (i_comment) {}
		// 
		// 	~ascii () {}
		// 	
		// 	std::string extension () {return ".dat";}
		// 	
		// 	virtual void open_file (std::string file_name, int file_type);
		// 
		// 	virtual void close_file ();
		// 
		// 	template <class datatype>
		// 	virtual void _write (std::string name, void *data) {
		// 		TRACE ("Writing...");
		// 		names.push_back (name);
		// 		double_data.resize (double_data.size () + 1);
		// 		float_data.resize (float_data.size () + 1);
		// 		int_data.resize (int_data.size () + 1);
		// 		if (type == &typeid (double)) {
		// 			double_data [(int) double_data.size () - 1].resize (n);
		// 			double_data [(int) double_data.size () - 1].assign ((double *) data, ((double *) data) + n);
		// 			types.push_back (&typeid (double));
		// 		} else if (type == &typeid (float)) {
		// 			float_data [(int) float_data.size () - 1].resize (n);
		// 			float_data [(int) float_data.size () - 1].assign ((float *) data, ((float *) data) + n);
		// 			types.push_back (&typeid (float));
		// 		} else if (type == &typeid (int)) {
		// 			int_data [(int) int_data.size () - 1].resize (n);
		// 			int_data [(int) int_data.size () - 1].assign ((int *) data, ((int *) data) + n);
		// 			types.push_back (&typeid (int));
		// 		}
		// 	}
		// 	virtual void write_scalar (std::string name, std::type_info *, void *) {
		// 		
		// 	}
		// 	virtual void read (std::string name, std::type_info *, void *) {
		// 		throw 0;
		// 	}
		// 	virtual void read_scalar (std::string name, std::type_info *, void *) {
		// 		throw 0;
		// 	}
		// 	
		// 	/*
		// 		TODO Write these methods
		// 	*/
		// 		
		// protected:
		// 	std::ofstream file_stream;
		// 	int n;
		// 	std::vector <const std::type_info *> types;
		// 	std::vector <std::vector <double>> double_data;
		// 	std::vector <std::vector <float>> float_data;
		// 	std::vector <std::vector <int>> int_data;
		// 	std::vector <std::string> names;
		// 	std::string comment;
		// };
		
		/*
			TODO Make legible headers
		*/
		
		// class netcdf : public format
		// {
		// public:
		// 	netcdf (int i_n, int i_n_offset = 0) :
		// 	n (i_n), n_offset (i_n_offset) {}
		// 
		// 	virtual ~netcdf () {}
		// 	
		// 	std::string extension () {return ".cdf";}
		// 	
		// 	virtual void open_file (std::string file_name, int file_type);
		// 
		// 	virtual void close_file ();
		// 
		// 	virtual void write (std::string name, double *);
		// 	virtual void write (std::string name, float *);
		// 	virtual void write (std::string name, int *);
		// 	
		// 	virtual void write_scalar (std::string name, double *);
		// 	virtual void write_scalar (std::string name, float *);
		// 	virtual void write_scalar (std::string name, int *);
		// 	
		// 	virtual void read (std::string name, double *);
		// 	virtual void read (std::string name, float *);
		// 	virtual void read (std::string name, int *);
		// 	
		// 	virtual void read_scalar (std::string name, double *);
		// 	virtual void read_scalar (std::string name, float *);
		// 	virtual void read_scalar (std::string name, int *);
		// 
		// protected:
		// 	int n;
		// 	int n_offset;
		// 	void * datafile;
		// 	netCDF::NcDim zDim;
		// };
	} /* one_d */
	
	namespace two_d
	{
		class virtual_format
		{
		public:
			virtual_format () {}
		
			~virtual_format () {}
		
			static std::string extension () {return "";}
			
			static void open_file (std::string file_name, int file_type, int n_max, int m_max, int l_max) {
				virtual_dumps [file_name];
			}
			
			static void close_file (std::string file_name) {}
			
			template <class datatype>
			static void write (std::string file_name, std::string name, datatype *data, int n = 1, int m = 1, int l = 1, int n_offset = 0, int m_offset = 0, int l_offset = 0) {
				virtual_dumps [file_name].add_var <datatype> (name, n, m);
				virtual_dumps [file_name].put <datatype> (name, (datatype *) data, n, m);
			}
		
			template <class datatype>
			static void write_scalar (std::string file_name, std::string name, datatype *data) {
				virtual_dumps [file_name].add_var <datatype> (name);
				virtual_dumps [file_name].put <datatype> (name, (datatype *) data);
			}
		
			template <class datatype>
			static void read (std::string file_name, std::string name, datatype *data, int n = 1, int m = 1, int l = 1, int n_offset = 0, int m_offset = 0, int l_offset = 0) {
				virtual_dumps [file_name].get <datatype> (name, (datatype *) data, n, m);
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j < m; ++j) {
						if (*(data + i * m + j) != *(data + i * m + j)) {
							FATAL ("NaN read in.");
							throw 0;
						}
					}
				}
			}
		
			template <class datatype>
			static void read_scalar (std::string file_name, std::string name, datatype *data) {
				virtual_dumps [file_name].get <datatype> (name, (datatype *) data);
			}
		};
			
		
		class netcdf
		{
		public:
			netcdf () {}
		
			virtual ~netcdf () {}
		
			static std::string extension () {return ".cdf";}
			
			static void open_file (std::string file_name, int file_type, int n_max, int m_max, int l_max);
		
			static void close_file (std::string file_name);
			
			template <class datatype>
			static void write (std::string file_name, std::string name, datatype *data, int n = 1, int m = 1, int l = 1, int n_offset = 0, int m_offset = 0, int l_offset = 0) {
				std::vector <size_t> offsets = {(size_t) n_offset, (size_t) m_offset};
				std::vector <size_t> sizes = {(size_t) n, (size_t) m};
				netCDF::NcVar ncdata = files [file_name]->addVar (name.c_str (), netcdf_type (&typeid (datatype)), dims [file_name]);
				ncdata.setFill (true, NULL);
				ncdata.putVar (offsets, sizes, data);
			}
		
			template <class datatype>
			static void write_scalar (std::string file_name, std::string name, datatype *data) {
				std::vector <netCDF::NcDim> scalar_dims;
				std::vector <size_t> scalar_offset;
				netCDF::NcVar ncdata = files [file_name]->addVar (name.c_str (), netcdf_type (&typeid (datatype)), scalar_dims);
				ncdata.putVar (scalar_offset, data);
			}
		
			template <class datatype>
			static void read (std::string file_name, std::string name, datatype *data, int n = 1, int m = 1, int l = 1, int n_offset = 0, int m_offset = 0, int l_offset = 0) {
				try {
					std::vector <size_t> offsets = {(size_t) n_offset, (size_t) m_offset};
					std::vector <size_t> sizes = {(size_t) n, (size_t) m};
					netCDF::NcVar ncdata = files [file_name]->getVar (name.c_str ());
					if (ncdata.isNull ()) {
						throw 0;
					}
					ncdata.getVar (offsets, sizes, data);
					for (int i = 0; i < n; ++i) {
						for (int j = 0; j < m; ++j) {
							if (*(data + i * m + j) != *(data + i * m + j)) {
								FATAL ("NaN read in.");
								throw 0;
							}
						}
					}
				} catch (netCDF::exceptions::NcBadName &e) {
					failures [file_name].push_back (name);
					WARN ("Variable " << name << " not found in file");
				} catch (int &e) {
					failures [file_name].push_back (name);
					WARN ("Variable " << name << " not found in file");
				}
			}
		
			template <class datatype>
			static void read_scalar (std::string file_name, std::string name, datatype *data) {
				try {
					DEBUG ("Reading scalar..." << name << " " << data);
					std::vector <size_t> scalar_offset;
					netCDF::NcVar ncdata = files [file_name]->getVar (name.c_str ());
					if (ncdata.isNull ()) {
						throw 0;
					}
					ncdata.getVar (scalar_offset, data);
					DEBUG ("Read " << *data);
				} catch (netCDF::exceptions::NcBadName &e) {
					failures [file_name].push_back (name);
					WARN ("Variable " << name << " not found in file");
				} catch (int &e) {
					failures [file_name].push_back (name);
					WARN ("Variable " << name << " not found in file");
				}
			}
			
		protected:
			static std::map <std::string, netCDF::NcFile *> files;
			static std::map <std::string, std::vector <netCDF::NcDim>> dims;
			static std::map <std::string, std::vector <std::string>> failures;
		};
		
		/*
			TODO Accept 1D outputs
		*/
	} /* two_d */
} /* io */

#endif /* end of include guard: IO_HPP_C1E9B6EF */
