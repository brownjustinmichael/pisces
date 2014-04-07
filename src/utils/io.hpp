/*!***********************************************************************
 * \file io.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <vector>
#include <array>
#include <map>
#include <string>
#include <iostream>
#include <memory>
#include <typeinfo>
#include "exceptions.hpp"
#include "../config.hpp"
#include <yaml-cpp/yaml.h>
#include <netcdfcpp.h>
#include "utils.hpp"

#ifndef IO_HPP_C1E9B6EF
#define IO_HPP_C1E9B6EF

namespace io
{	
	class virtual_dump
	{
	public:
		virtual_dump () {}
		
		virtual ~virtual_dump () {}
		
		template <class datatype>
		typename std::map <std::string, std::vector <datatype>>::iterator begin ();
		
		template <class datatype>
		typename std::map <std::string, std::vector <datatype>>::iterator end ();
		
		template <class datatype>
		void add_var (std::string name, int n = 1, int m = 1) {
			TRACE ("ADDING VAR");
			if (typeid (datatype) == typeid (double)) {
				double_map [name].resize (n * m);
			} else if (typeid (datatype) == typeid (float)) {
				float_map [name].resize (n * m);
			} else if (typeid (datatype) == typeid (int)) {
				int_map [name].resize (n * m);
			} else {
				ERROR ("Unknown type");
				throw 0;
			}
			dims [name] [0] = n;
			dims [name] [1] = m;
			TRACE ("Done.");
		}
		
		template <class datatype>
		void put (std::string name, void *data, int n = 1, int m = 1, int ldm = -1) {
			TRACE ("Putting " << name);
			ldm = ldm == -1 ? m : ldm;
			if (typeid (datatype) == typeid (double)) {
				utils::matrix_copy (m, n, (double *) data, &(double_map [name] [0]), ldm, dims [name] [1]);
			} else if (typeid (datatype) == typeid (float)) {
				utils::matrix_copy (m, n, (float *) data, &(float_map [name] [0]), ldm, dims [name] [1]);
			} else if (typeid (datatype) == typeid (int)) {
				utils::matrix_copy (m, n, (int *) data, &(int_map [name] [0]), ldm, dims [name] [1]);
			} else {
				ERROR ("Unknown type");
				throw 0;
			}
			TRACE ("Done.");
		}
		
		template <class datatype>
		void get (std::string name, void *data, int n = 1, int m = 1, int ldm = -1) {
			TRACE ("Getting...");
			ldm = ldm == -1 ? m : ldm;
			if (typeid (datatype) == typeid (double)) {
				utils::matrix_copy (m, n, &(double_map [name] [0]), (double *) data, dims [name] [1], ldm);
			} else if (typeid (datatype) == typeid (float)) {
				utils::matrix_copy (m, n, &(float_map [name] [0]), (float *) data, dims [name] [1], ldm);
			} else if (typeid (datatype) == typeid (int)) {
				utils::matrix_copy (m, n, &(int_map [name] [0]), (int *) data, dims [name] [1], ldm);
			} else {
				ERROR ("Unknown type");
				throw 0;
			}
			TRACE ("Done.");
		}
		
		template <class datatype>
		datatype &index (std::string name, int i = 0, int j = 0);
			
		std::map <std::string, std::array <int, 2>> dims;
	private:
		std::map <std::string, std::vector <int>> int_map;
		std::map <std::string, std::vector <float>> float_map;
		std::map <std::string, std::vector <double>> double_map;
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
	
	
	class format
	{
	public:
		virtual ~format () {}
		
		virtual std::string extension () = 0;
		
		virtual void to_file (const char *file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs) = 0;
		
		virtual void from_file (const char *file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs) = 0;
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
		}
		
		template <class datatype>
		void append_scalar (std::string name, datatype *data_ptr) {
			TRACE ("Appending " << name << " to output...");
			for (int i = 0; i < (int) scalar_names.size (); ++i) {
				if (scalar_names [i] == name) {
					WARN ("Reuse of name " << name);
					scalar_types [i] = &typeid (datatype);
					scalar_ptrs [i] = (void *) data_ptr;
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
		 * This function should be overwritten by subclasses, though it may 
		 * contain a call to this function, which will output with default 
		 * datatype representation in C++.
		 *********************************************************************/
		virtual void to_file () {
			TRACE ("Sending to file...");
			format_ptr->to_file (file_name.c_str (), (int) names.size (), &names [0], &types [0], &data_ptrs [0], (int) scalar_names.size (), &scalar_names [0], &scalar_types [0], &scalar_ptrs [0]);
		}
		
	protected:
		std::shared_ptr <format> format_ptr;
		std::string file_name;
		std::vector <std::string> names;
		std::vector <std::string> scalar_names;
		std::vector <const std::type_info*> types;
		std::vector <const std::type_info*> scalar_types;
		std::vector <void *> data_ptrs; //!< A vector of integer pointers to the arrays of data
		std::vector <void *> scalar_ptrs; //!< A vector of integer pointers to the arrays of data
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
			// TRACE ("Gathering from file...");
			format_ptr->from_file (file_name.c_str (), (int) names.size (), &names [0], &types [0], &data_ptrs [0], (int) scalar_names.size (), &scalar_names [0], &scalar_types [0], &scalar_ptrs [0]);
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
			ascii (int i_n) : 
			n (i_n) {}
		
			~ascii () {}
			
			std::string extension () {return ".dat";}
		
			/*!*******************************************************************
			 * \brief Outputs to file_name
			 *********************************************************************/
			void to_file (const char *file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs);
		
			void from_file (const char *file_name, int n_data_ptrs = 0, std::string *names = NULL, const std::type_info **types = NULL, void **data_ptrs = NULL, int n_scalar_ptrs = 0, std::string *scalar_names = NULL, const std::type_info **scalar_types = NULL, void **scalar_ptrs = NULL);
		
		protected:
			int n;
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
			
			void to_file (const char *file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs);
	
			void from_file (const char *file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs);
	
		protected:
			int n;
		};
	} /* one_d */
	
	namespace two_d
	{
		class virtual_format : public format
		{
		public:
			virtual_format (virtual_dump *i_dump, int i_n, int i_m) :
			dump (i_dump), n (i_n), m (i_m) {}
		
			~virtual_format () {}
		
			std::string extension () {return "";}
		
			void to_file (const char *file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs) {
				std::vector <std::string> failures;
				for (int i = 0; i < n_data_ptrs; ++i) {
					if (*types [i] == typeid (double)) {
						dump->add_var <double> (names [i], n, m);
						dump->put <double> (names [i], (double *) data_ptrs [i], n, m);
					} else if (*types [i] == typeid (float)) {
						dump->add_var <float> (names [i], n, m);
						dump->put <float> (names [i], (float *) data_ptrs [i], n, m);
					} else if (*types [i] == typeid (int)) {
						dump->add_var <int> (names [i], n, m);
						dump->put <int> (names [i], (int *) data_ptrs [i], n, m);
					} else {
						throw 0;
					}
				}
				
				for (int i = 0; i < n_scalar_ptrs; ++i) {
					if (*scalar_types [i] == typeid (double)) {
						dump->add_var <double> (scalar_names [i]);
						dump->put <double> (scalar_names [i], (double *) scalar_ptrs [i]);
					} else if (*scalar_types [i] == typeid (float)) {
						dump->add_var <float> (scalar_names [i]);
						dump->put <float> (scalar_names [i], (float *) scalar_ptrs [i]);
					} else if (*scalar_types [i] == typeid (int)) {
						dump->add_var <int> (scalar_names [i]);
						dump->put <int> (scalar_names [i], (int *) scalar_ptrs [i]);
					} else {
						throw 0;
					}
				}
			}
			
			void from_file (const char *file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs) {
				for (int i = 0; i < n_data_ptrs; ++i) {
					if (*types [i] == typeid (double)) {
						dump->get <double> (names [i], (double *) data_ptrs [i], n, m);
					} else if (*types [i] == typeid (float)) {
						dump->get <float> (names [i], (float *) data_ptrs [i], n, m);
					} else if (*types [i] == typeid (int)) {
						dump->get <int> (names [i], (int *) data_ptrs [i], n, m);
					} else {
						throw 0;
					}
				}
				
				for (int i = 0; i < n_scalar_ptrs; ++i) {
					if (*scalar_types [i] == typeid (double)) {
						dump->get <double> (scalar_names [i], (double *) scalar_ptrs [i]);
					} else if (*scalar_types [i] == typeid (float)) {
						dump->get <float> (scalar_names [i], (float *) scalar_ptrs [i]);
					} else if (*scalar_types [i] == typeid (int)) {
						dump->get <int> (scalar_names [i], (int *) scalar_ptrs [i]);
					} else {
						throw 0;
					}
				}
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
		
			void to_file (const char *file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs);

			void from_file (const char *file_name, int n_data_ptrs = 0, std::string *names = NULL, const std::type_info **types = NULL, void **data_ptrs = NULL, int n_scalar_ptrs = 0, std::string *scalar_names = NULL, const std::type_info **scalar_types = NULL, void **scalar_ptrs = NULL);
	
		protected:
			int n, m;
			std::vector <std::string> failures;
		};
		
		/*
			TODO Accept 1D outputs
		*/
	} /* two_d */
} /* io */

#endif /* end of include guard: IO_HPP_C1E9B6EF */
