/*!***********************************************************************
 * \file io.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "io.hpp"
#include "exceptions.hpp"
#include "../config.hpp"

namespace io
{
	simple_header::simple_header (int i_n_data_headers, std::string *i_data_headers, std::string i_comment_string) : header () {
		int i;

		TRACE ("Constructing...");

		header_string = i_comment_string;
			
		for (i = 0; i < i_n_data_headers; ++i)
		{
			header_string.append (" ");
			header_string.append (i_data_headers [i]);
		}
		
		header_string.append ("\n");

		TRACE ("Constructed.");
	}
	
	output::output (header *i_header_ptr, int i_n, int i_n_data_ptrs, double **i_data_ptrs) {
		int i;
		n = i_n;
		n_data_ptrs = i_n_data_ptrs;
		header_ptr = i_header_ptr;
		
		TRACE ("Constructing...");

		for (i = 0; i < i_n_data_ptrs; ++i)
		{
			double_ptrs.push_back (i_data_ptrs [i]);
			int_ptrs.push_back (NULL);
		}

		TRACE ("Constructed.");
	}
	
	void output::append (double *data_ptr) {
		TRACE ("Appending double to output...");
		
		++n_data_ptrs;
		double_ptrs.push_back (data_ptr);
		int_ptrs.push_back (NULL);
		
		TRACE ("Appended.");
	}
	
	void output::append (int *data_ptr) {
		TRACE ("Appending int to output...");
		
		++n_data_ptrs;
		double_ptrs.push_back (NULL);
		int_ptrs.push_back (data_ptr);
		
		TRACE ("Appended.");
	}
	
	void output::simple_to_file (std::string file_name) {
		int i, j;
		std::ofstream file_stream; // A file stream object to be used when writing to file 
		
		TRACE ("Beginning output...");		
		
		INFO ("Outputting to file " << file_name << "...");
		
		file_stream.open (file_name.c_str ());
		
		if (! file_stream.is_open ()) {
			exceptions::file_exception failure;
			ERROR ("Failed to open file " << file_name);
			throw failure;
		}
		
		file_stream << header_ptr->output_header ();
						
		for (i = 0; i < n; ++i)
		{
			for (j = 0; j < n_data_ptrs; ++j)
			{
				if (double_ptrs [j] == NULL) {
					file_stream << int_ptrs [j] [i];
				} else {
					file_stream << double_ptrs [j] [i];
				}
				file_stream << ' ';
			}
			file_stream << '\n';
		}
				
		file_stream.close ();
		
		TRACE ("Output to file.");
	}
	
	std::string incremental_output::generate_file_name () {

		TRACE ("Generating new file_name...");

		if (n_outputs >= std::pow (10., int_width)) {
			ERROR ("Exceeded maximum number of outputs")
		}
		std::ostringstream new_file_name;
		new_file_name << file_base << std::setfill ('0') << std::setw (int_width) << n_outputs << file_extension;

		++n_outputs;

		TRACE ("New file name is " << new_file_name.str () << ".");
		
		return std::string (new_file_name.str ());
	}
} /* io */