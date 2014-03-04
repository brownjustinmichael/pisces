/*!***********************************************************************
 * \file io.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <memory>
#include <typeinfo>
#include "exceptions.hpp"
#include "../config.hpp"
#include <yaml-cpp/yaml.h>


#ifndef IO_HPP_C1E9B6EF
#define IO_HPP_C1E9B6EF

namespace io
{	
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
	
	
	class format
	{
	public:
		format () {}
		
		virtual ~format () {
			// printf ("Destroying format\n");
		}
		
		virtual std::string extension () = 0;
		
		virtual void to_file (std::string file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs) = 0;
		
		virtual void from_file (std::string file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs) = 0;
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
		file_name (i_file_name) {};
		
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
			types.push_back (&typeid (datatype));
			names.push_back (name);
			data_ptrs.push_back ((void *) data_ptr);
		}
		
		template <class datatype>
		void append_scalar (std::string name, datatype *data_ptr) {
			TRACE ("Appending " << name << " to output...");
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
			format_ptr->to_file (file_name + format_ptr->extension (), (int) names.size (), &names [0], &types [0], &data_ptrs [0], (int) scalar_names.size (), &scalar_names [0], &scalar_types [0], &scalar_ptrs [0]);
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
		file_format (i_file_format),
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
		file_name (i_file_name) {};
		
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
			TRACE ("Appending " << name << " to input...");
			types.push_back (&typeid (datatype));
			names.push_back (name);
			data_ptrs.push_back ((void *) data_ptr);
		}
		
		template <class datatype>
		void append_scalar (std::string name, datatype *data_ptr) {
			TRACE ("Appending " << name << " to output...");
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
			TRACE ("Sending to file...");
			format_ptr->from_file (file_name + format_ptr->extension (), (int) names.size (), &names [0], &types [0], &data_ptrs [0], (int) scalar_names.size (), &scalar_names [0], &scalar_types [0], &scalar_ptrs [0]);
		}
		
	protected:
		std::shared_ptr <format> format_ptr;
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
			void to_file (std::string file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs);
		
			void from_file (std::string file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs);
		
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
			
			void to_file (std::string file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs);

			void from_file (std::string file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs);
	
		protected:
			int n;
		};
	} /* one_d */
	
	namespace two_d
	{
		class netcdf : public format
		{
		public:
			netcdf (int i_n = 0, int i_m = 0) :
			n (i_n), m (i_m) {
				// if (n == 0)
			}
		
			virtual ~netcdf () {}
		
			std::string extension () {return ".cdf";}
		
			void to_file (std::string file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs);

			void from_file (std::string file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs);
	
		protected:
			int n, m;
		};
		
		/*
			TODO Accept 1D outputs
		*/
	} /* two_d */
} /* io */

#endif /* end of include guard: IO_HPP_C1E9B6EF */
