// 
//! \file io.cpp
//  io
//  
//  Created by Justin Brown on 2013-03-26.
//  Copyright 2013 Justin Brown. All rights reserved.
// 

#include <cmath>
#include <iomanip>
#include "io.hpp"
#include "exceptions.hpp"
#include "../config.hpp"

namespace io
{
	simple_header::simple_header (int i_n_data_headers, std::string *i_data_headers, std::string i_comment_string) : header () {
		int i;

		LOG4CXX_TRACE (config::logger, "Constructing...");

		header_string = i_comment_string;
			
		for (i = 0; i < i_n_data_headers; ++i)
		{
			header_string.append (" ");
			header_string.append (i_data_headers [i]);
		}
		
		header_string.append ("\n");

		LOG4CXX_TRACE (config::logger, "Constructed.");
	}
	
	output::output (header *i_header_ptr, int i_n, int i_n_data_ptrs, double **i_data_ptrs) {
		int i;
		n = i_n;
		n_data_ptrs = i_n_data_ptrs;
		header_ptr = i_header_ptr;
		
		LOG4CXX_TRACE (config::logger, "Constructing...");

		for (i = 0; i < i_n_data_ptrs; ++i)
		{
			double_ptrs.push_back (i_data_ptrs [i]);
			int_ptrs.push_back (NULL);
		}

		LOG4CXX_TRACE (config::logger, "Constructed.");
	}
	
	void output::append (double *data_ptr) {
		++n_data_ptrs;
		double_ptrs.push_back (data_ptr);
		int_ptrs.push_back (NULL);
	}
	
	void output::append (int *data_ptr) {
		++n_data_ptrs;
		double_ptrs.push_back (NULL);
		int_ptrs.push_back (data_ptr);
	}
	
	void output::to_file (std::string file_name) {
		int i, j;
		
		LOG4CXX_TRACE (config::logger, "Beginning output...");		
		
		LOG4CXX_INFO (config::logger, "Outputting to file " << file_name << "...");
		
		file_stream.open (file_name.c_str ());
		
		if (not file_stream.is_open ()) {
			exceptions::file_exception failure;
			LOG4CXX_ERROR (config::logger, "Failed to open file " << file_name);
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
		
		LOG4CXX_TRACE (config::logger, "Output to file.");
	}
	
	std::string incremental_output::generate_file_name () {

		LOG4CXX_TRACE (config::logger, "Generating new file_name...");

		if (n_outputs >= std::pow (10, int_width)) {
			LOG4CXX_ERROR (config::logger, "Exceeded maximum number of outputs")
		}
		std::ostringstream new_file_name;
		new_file_name << file_base << std::setfill ('0') << std::setw (int_width) << n_outputs << file_extension;

		++n_outputs;

		LOG4CXX_TRACE (config::logger, "New file name is " << new_file_name.str () << ".");
		return std::string (new_file_name.str ());
	}
} /* io */