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
#include "../config.hpp"

#ifndef IO_HPP_C1E9B6EF
#define IO_HPP_C1E9B6EF

namespace io
{
	template <class datatype>
	class parameters
	{
	public:
		std::string filename;

		datatype diffusion_coeff;
		datatype advection_coeff;
		int timesteps;
		int gridpoints;
		int output_every;
		int n_iterations;
		int n, nmp, nrhs, nb;
		datatype max_timestep, courant_factor;
		datatype scale, width, mean, sigma;

		parameters (std::string i_filename);
		
		virtual ~parameters () {}
	};
	
	class format
	{
	public:
		format () {}
		
		virtual ~format () {
			printf ("Destroying format\n");
		}
		
		virtual void to_file (std::string file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs) = 0;
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
			printf ("Destroying output\n");
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
		
		/*!*******************************************************************
		 * \brief A virtual function to output the data to file
		 * 
		 * This function should be overwritten by subclasses, though it may 
		 * contain a call to this function, which will output with default 
		 * datatype representation in C++.
		 *********************************************************************/
		virtual void to_file () {
			format_ptr->to_file (file_name, (int) names.size (), &names [0], &types [0], &data_ptrs [0]);
		}
		
	protected:
		std::shared_ptr <format> format_ptr;
		std::string file_name;
		std::vector <std::string> names;
		std::vector <const std::type_info*> types;
		std::vector <void *> data_ptrs; //!< A vector of integer pointers to the arrays of data
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
			if (count % output_every == 0) {
				char buffer [100];
				snprintf (buffer, 100, file_format.c_str (), count / output_every);
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
		
			/*!*******************************************************************
			 * \brief Outputs to file_name
			 *********************************************************************/
			void to_file (std::string file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs);
		
		protected:
			int n;
		};
		
		/*
			TODO Make legible headers
		*/
	} /* one_d */
	
	namespace two_d
	{
		class netcdf : public format
		{
		public:
			netcdf (int i_n, int i_m) :
			n (i_n), m (i_m) {}
		
			virtual ~netcdf () {}
		
			void to_file (std::string file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs);
	
		protected:
			int n, m;
		};
		
		/*
			TODO Accept 1D outputs
		*/
	} /* two_d */
} /* io */

#endif /* end of include guard: IO_HPP_C1E9B6EF */
