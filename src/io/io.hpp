// 
//! \file io.hpp
//  io
//  
//  Created by Justin Brown on 2013-03-25.
//  Copyright 2013 Justin Brown. All rights reserved.
// 

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "../config.hpp"

#ifndef IO_HPP_C1E9B6EF
#define IO_HPP_C1E9B6EF

namespace io
{
	//! \brief A base class that is essentially a null header
	//
	//! Subclasses inheriting from this class should overload the output_header method, which outputs a string that will represent the header of the output file.
	class header
	{
	public:
		//! \brief Method that generates the header for an output_stream
		//
		//! This method should be overloaded by subclasses to be replaced with more elaborate headers.
		//! \return A string which will be used as the header in the associated output_stream
		virtual std::string output_header () {return std::string ("");}
		// virtual ~header ();
	};
	
	//! \brief A simpler header class
	//
	//! This class generates a simple header that contains just the names of the columns
	class simple_header : public header
	{
	private:
		std::string header_string; //!< The string form of the header
	public:
		//! \param i_n_data_headers An integer number of column headers to be read from i_data_headers
		//! \param i_data_headers An array of strings that are the column headers
		//! \param i_comment_char A string that will be used at the beginning of each commented line
		simple_header (int i_n_data_headers, std::string *i_data_headers, std::string i_comment_string = "#");
		using header::output_header;
		std::string output_header () {return header_string;}
	};
	
	// \brief An abstract output stream base class that generates output files
	class output_stream
	{
	protected:
		int n; //!< An integer number of elements in each array (regardless of dimension)
		int n_data_ptrs; //!< An integer number of arrays to output
		std::vector<double *> data_ptrs; //!< A vector of double pointers to the arrays of data
		header *header_ptr; //!< A pointer to a header object, which contains the details regarding the construction of the header
		std::ofstream file_stream; //!< A file stream object to be used when writing to file 
	public:
		output_stream (header *i_header_ptr, int i_n, int i_n_data_ptrs, double **i_data_ptrs);
		virtual void output (std::string i_file_name);
		// virtual ~output_stream ();
	};
	
	//! \brief A simple 1D implementation of the output_stream class
	//
	//! This class is a simple implementation of the output_stream class that can only deal with 1D data
	class simple_output_stream_1D : public output_stream
	{
	private:
		std::string file_name; //!< A string containing the file name where the class should output
	public:
		simple_output_stream_1D (std::string i_file_name, int i_n, int i_n_data_ptrs, double **i_data_ptrs) : output_stream (new header, i_n, i_n_data_ptrs, i_data_ptrs) {file_name = i_file_name;}
		using output_stream::output;
		void output () {output_stream::output (file_name);} 
	};
	
	//! \brief An incremental 1D implementation of the output_stream class
	//
	//! This class is an implementation of the output_stream class that can only deal with 1D data but that increments the file number with each iteration
	class incremental_output_stream_1D : public output_stream
	{
	private:
		int n_outputs; //!< An incremental integer that will be appended to the end of each file_base
		int int_width; //!< The total number of characters the integer in the file name can have
		std::string file_base;
		std::string file_extension;
	public:
		incremental_output_stream_1D (std::string i_file_base, std::string i_file_extension, int i_int_width, header *i_header_ptr, int i_n, int i_n_data_ptrs, double **i_data_ptrs) : output_stream (i_header_ptr, i_n, i_n_data_ptrs, i_data_ptrs) {n_outputs = 0; int_width = i_int_width, file_base = i_file_base; file_extension = i_file_extension;}
		std::string generate_file_name ();
		void output () {output_stream::output (generate_file_name ());} 
	};
} /* io */

#endif /* end of include guard: IO_HPP_C1E9B6EF */
