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
#include <netcdfcpp.h>
#include "utils.hpp"

namespace io
{	
	class virtual_dump
	{
	public:
		virtual_dump () {}
		
		virtual ~virtual_dump () {
			for (std::map <std::string, void *>::iterator iter = data_map.begin (); iter != data_map.end (); iter++) {
				free (iter->second);
			}
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
		parameters (std::string filename) {
			YAML::Node::operator= (YAML::LoadFile (filename));
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
	
	class format
	{
	public:
		virtual ~format () {}
		
		virtual std::string extension () = 0;
		
		virtual void open_file (std::string file_name, int file_type) = 0;
		
		virtual void close_file () = 0;
		
		virtual void write (std::string name, double *) = 0;
		virtual void write (std::string name, float *) = 0;
		virtual void write (std::string name, int *) = 0;
		
		virtual void write_scalar (std::string name, double *) = 0;
		virtual void write_scalar (std::string name, float *) = 0;
		virtual void write_scalar (std::string name, int *) = 0;
		
		virtual void read (std::string name, double *) = 0;
		virtual void read (std::string name, float *) = 0;
		virtual void read (std::string name, int *) = 0;
		
		virtual void read_scalar (std::string name, double *) = 0;
		virtual void read_scalar (std::string name, float *) = 0;
		virtual void read_scalar (std::string name, int *) = 0;
		
		const static int read_file = 0;
		const static int replace_file = 1;
		const static int append_file = 2;
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
		output (format* i_format_ptr, std::string i_file_name = "out") :
		format_ptr (i_format_ptr),
		file_name (i_file_name + format_ptr->extension ()) {};
		
		virtual ~output () {
			// printf ("Destroying output\n");
		}
		
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
			TRACE ("Scalar appended.");
		}
		
		/*!*******************************************************************
		 * \brief A virtual function to output the data to file
		 * 
		 * This function should be overwritten by subclasses, though it may 
		 * contain a call to this function, which will output with default 
		 * datatype representation in C++.
		 *********************************************************************/
		virtual void to_file () {
			TRACE ("Sending to file...");
			
			INFO ("Outputting to file " << file_name << "...");
			
			// DEBUG ("Scalars");
			format_ptr->open_file (file_name.c_str (), format::replace_file);
			for (int i = 0; i < (int) scalar_names.size (); ++i) {
				if (scalar_types [i] == &typeid (double)) {
					format_ptr->write_scalar (scalar_names [i], (double *) scalar_ptrs [i]);
				} else if (scalar_types [i] == &typeid (float)) {
					format_ptr->write_scalar (scalar_names [i], (float *) scalar_ptrs [i]);
				} else if (scalar_types [i] == &typeid (int)) {
					format_ptr->write_scalar (scalar_names [i], (int *) scalar_ptrs [i]);
				}
			}

			// DEBUG ("Data");
			for (int i = 0; i < (int) names.size (); ++i) {
				if (types [i] == &typeid (double)) {
					format_ptr->write (names [i], (double *) data_ptrs [i]);
				} else if (types [i] == &typeid (float)) {
					format_ptr->write (names [i], (float *) data_ptrs [i]);
				} else if (types [i] == &typeid (int)) {
					format_ptr->write (names [i], (int *) data_ptrs [i]);
				}
			}
			
			// DEBUG ("Functions");
			for (int i = 0; i < (int) functor_names.size (); ++i) {
				if (functor_types [i] == &typeid (double)) {
					format_ptr->write (functor_names [i], ((format_functor <double> *) functor_ptrs [i])->calculate ());
				} else if (functor_types [i] == &typeid (float)) {
					format_ptr->write (functor_names [i], ((format_functor <float> *) functor_ptrs [i])->calculate ());
				} else if (functor_types [i] == &typeid (int)) {
					format_ptr->write (functor_names [i], ((format_functor <int> *) functor_ptrs [i])->calculate ());
				}
			}
			format_ptr->close_file ();
		}
		
	protected:
		std::shared_ptr <format> format_ptr;
		std::string file_name;
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
	
	class incremental : public output
	{
	public:
		incremental (format* i_format_ptr, std::string i_file_format, int i_output_every = 1) :
		output (i_format_ptr),
		file_format (i_file_format + format_ptr->extension ()),
		output_every (i_output_every > 0 ? i_output_every : 1),
		count (0) {}
		
		virtual ~incremental () {}
	
		void to_file () {
			TRACE ("Sending to file...");
			if (count % output_every == 0) {
				char buffer [file_format.size () * 2];
				snprintf (buffer, file_format.size () * 2, file_format.c_str (), count / output_every);
				file_name = buffer;
				output::to_file ();
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
		input (format* i_format_ptr, std::string i_file_name = "in") :
		format_ptr (i_format_ptr),
		file_name (i_file_name + format_ptr->extension ()) {};
		
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
		virtual void from_file () {
			INFO ("Inputting from file " << file_name << "...");

			format_ptr->open_file (file_name.c_str (), format::read_file);
			for (int i = 0; i < (int) names.size (); ++i) {
				if (types [i] == &typeid (double)) {
					format_ptr->read (names [i], (double *) data_ptrs [i]);
				} else if (types [i] == &typeid (float)) {
					format_ptr->read (names [i], (float *) data_ptrs [i]);
				} else if (types [i] == &typeid (int)) {
					format_ptr->read (names [i], (int *) data_ptrs [i]);
				}
			}
			
			for (int i = 0; i < (int) scalar_names.size (); ++i) {
				if (scalar_types [i] == &typeid (double)) {
					format_ptr->read_scalar (scalar_names [i], (double *) scalar_ptrs [i]);
				} else if (scalar_types [i] == &typeid (float)) {
					format_ptr->read_scalar (scalar_names [i], (float *) scalar_ptrs [i]);
				} else if (scalar_types [i] == &typeid (int)) {
					format_ptr->read_scalar (scalar_names [i], (int *) scalar_ptrs [i]);
				}
			}
			format_ptr->close_file ();
		}
		
	protected:
		format *format_ptr;
		std::string file_name;
		std::vector <std::string> names;
		std::vector <const std::type_info*> types;
		std::vector <void *> data_ptrs; //!< A vector of integer pointers to the arrays of data
		std::vector <std::string> scalar_names;
		std::vector <const std::type_info*> scalar_types;
		std::vector <void *> scalar_ptrs; //!< A vector of integer pointers to the arrays of data
	};
	
	namespace one_d
	{
		/*!*******************************************************************
		 * \brief A simple implementation of the output class
		 * 
		 * This class is a simple implementation of the output class.
		 *********************************************************************/
		class ascii : public format
		{
		public:
			/*!*******************************************************************
			 * \param i_file_name The string name of file for output
			 * \param i_n The integer number of points in the data
			 * \param i_output_every An integer number of steps between outputs
			 *********************************************************************/
			ascii (int i_n, std::string i_comment = "#") : 
			n (i_n),
			comment (i_comment) {}
		
			~ascii () {}
			
			std::string extension () {return ".dat";}
			
			virtual void open_file (std::string file_name, int file_type);
		
			virtual void close_file ();
		
			virtual void write (std::string name, double *);
			virtual void write (std::string name, float *);
			virtual void write (std::string name, int *);
			
			virtual void write_scalar (std::string name, double *);
			virtual void write_scalar (std::string name, float *);
			virtual void write_scalar (std::string name, int *);
			
			virtual void read (std::string name, double *);
			virtual void read (std::string name, float *);
			virtual void read (std::string name, int *);
			
			virtual void read_scalar (std::string name, double *);
			virtual void read_scalar (std::string name, float *);
			virtual void read_scalar (std::string name, int *);
			
			/*
				TODO Write these methods
			*/
				
		protected:
			std::ofstream file_stream;
			int n;
			std::vector <const std::type_info *> types;
			std::vector <std::vector <double>> double_data;
			std::vector <std::vector <float>> float_data;
			std::vector <std::vector <int>> int_data;
			std::vector <std::string> names;
			std::string comment;
		};
		
		/*
			TODO Make legible headers
		*/
		
		class netcdf : public format
		{
		public:
			netcdf (int i_n) :
			n (i_n) {}
		
			virtual ~netcdf () {}
			
			std::string extension () {return ".cdf";}
			
			virtual void open_file (std::string file_name, int file_type);
		
			virtual void close_file ();
		
			virtual void write (std::string name, double *);
			virtual void write (std::string name, float *);
			virtual void write (std::string name, int *);
			
			virtual void write_scalar (std::string name, double *);
			virtual void write_scalar (std::string name, float *);
			virtual void write_scalar (std::string name, int *);
			
			virtual void read (std::string name, double *);
			virtual void read (std::string name, float *);
			virtual void read (std::string name, int *);
			
			virtual void read_scalar (std::string name, double *);
			virtual void read_scalar (std::string name, float *);
			virtual void read_scalar (std::string name, int *);
		
		protected:
			int n;
			void * datafile, *zDim;
		};
	} /* one_d */
	
	namespace two_d
	{
		class average_profile : public format
		{
		public:
			average_profile (format *i_format, int i_n, int i_m) :
			base_format (i_format),
			n (i_n),
			m (i_m) {
				double_buffer.resize (m);
				int_buffer.resize (m);
				float_buffer.resize (m);
			}
			
			virtual ~average_profile () {}
			
			std::string extension () {return base_format->extension ();}
			
			void open_file (std::string file_name, int file_type) {
				TRACE ("Opening file...");
				base_format->open_file (file_name, file_type);
			}
			
			void close_file () {
				base_format->close_file ();
			}
			
			void write (std::string name, double *data) {
				TRACE ("Making average...");
				for (int j = 0; j < m; ++j) {
					double_buffer [j] = 0.0;
					for (int i = 0; i < n; ++i) {
						double_buffer [j] += data [i * m + j];
					}
					double_buffer [j] /= n;
				}
				base_format->write (name, &double_buffer [0]);
			}
			void write (std::string name, float *data) {
				TRACE ("Making average...");
				for (int j = 0; j < m; ++j) {
					float_buffer [j] = 0.0;
					for (int i = 0; i < n; ++i) {
						float_buffer [j] += data [i * m + j];
					}
					float_buffer [j] /= n;
				}
				base_format->write (name, &float_buffer [0]);
			}
			void write (std::string name, int *data) {
				TRACE ("Making average...");
				for (int j = 0; j < m; ++j) {
					int_buffer [j] = 0;
					for (int i = 0; i < n; ++i) {
						int_buffer [j] += data [i * m + j];
					}
					int_buffer [j] /= n;
				}
				base_format->write (name, &int_buffer [0]);
			}
			
			void write_scalar (std::string name, double *data) {
				base_format->write_scalar (name, data);
			}
			void write_scalar (std::string name, float *data) {
				base_format->write_scalar (name, data);
			}
			void write_scalar (std::string name, int *data) {
				base_format->write_scalar (name, data);
			}
			
			void read (std::string name, double *data) {
				throw 0;
			}
			void read (std::string name, float *data) {
				throw 0;
			}
			void read (std::string name, int *data) {
				throw 0;
			}
			
			void read_scalar (std::string name, double *data) {
				throw 0;
			}
			void read_scalar (std::string name, float *data) {
				throw 0;
			}
			void read_scalar (std::string name, int *data) {
				throw 0;
			}
			
		private:
			format *base_format;
			int n, m;
			std::vector <int> int_buffer;
			std::vector <double> double_buffer;
			std::vector <float> float_buffer;
		};
		
		class virtual_format : public format
		{
		public:
			virtual_format (virtual_dump *i_dump, int i_n, int i_m) :
			dump (i_dump), n (i_n), m (i_m) {}
		
			~virtual_format () {}
		
			std::string extension () {return "";}
			
			void open_file (std::string name, int file_type) {}
			void close_file () {}
			
			void write (std::string name, double *data) {
				dump->add_var <double> (name, n, m);
				dump->put <double> (name, (double *) data, n, m);
			}
			void write (std::string name, float *data) {
				dump->add_var <float> (name, n, m);
				dump->put <float> (name, (float *) data, n, m);
			}
			void write (std::string name, int *data) {
				dump->add_var <int> (name, n, m);
				dump->put <int> (name, (int *) data, n, m);
			}
			
			void write_scalar (std::string name, double *data) {
				dump->add_var <double> (name);
				dump->put <double> (name, (double *) data);
			}
			void write_scalar (std::string name, float *data) {
				dump->add_var <float> (name);
				dump->put <float> (name, (float *) data);
			}
			void write_scalar (std::string name, int *data) {
				dump->add_var <int> (name);
				dump->put <int> (name, (int *) data);
			}
			
			void read (std::string name, double *data) {
				dump->get <double> (name, (double *) data, n, m);
			}
			void read (std::string name, float *data) {
				dump->get <float> (name, (float *) data, n, m);
			}
			void read (std::string name, int *data) {
				dump->get <int> (name, (int *) data, n, m);
			}
			
			void read_scalar (std::string name, double *data) {
				dump->get <double> (name, (double *) data);
			}
			void read_scalar (std::string name, float *data) {
				dump->get <float> (name, (float *) data);
			}
			void read_scalar (std::string name, int *data) {
				dump->get <int> (name, (int *) data);
			}
		
		private:
			virtual_dump *dump;
			int n, m;
		};
			
		
		class netcdf : public format
		{
		public:
			netcdf (int i_n = 0, int i_m = 0) :
			n (i_n), m (i_m) {
				// if (n == 0)
			}
		
			virtual ~netcdf () {}
		
			std::string extension () {return ".cdf";}
			
			virtual void open_file (std::string file_name, int file_type);
		
			virtual void close_file ();
			
			void write (std::string name, double *data);
			void write (std::string name, float *data);
			void write (std::string name, int *data);
			
			void write_scalar (std::string name, double *data);
			void write_scalar (std::string name, float *data);
			void write_scalar (std::string name, int *data);
			
			void read (std::string name, double *data);
			void read (std::string name, float *data);
			void read (std::string name, int *data);
			
			void read_scalar (std::string name, double *data);
			void read_scalar (std::string name, float *data);
			void read_scalar (std::string name, int *data);
			
		protected:
			int n, m;
			std::vector <std::string> failures;
			void *datafile, *zDim, *xDim;
		};
		
		/*
			TODO Accept 1D outputs
		*/
	} /* two_d */
} /* io */

#endif /* end of include guard: IO_HPP_C1E9B6EF */
