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
	/****************** BEGIN: read_parameter *******************/

	read_params_txt::read_params_txt (std::string i_filename)
	{
		std::string temp;
		filename = i_filename;

		std::ifstream input_stream (filename);

		if (input_stream.is_open())
		{
			input_stream >> temp;
			input_stream >> inputParam[temp].asDouble;		// diffusion_coeff 
			input_stream >> temp;
			input_stream >> inputParam[temp].asDouble;		// advection_coeff
			input_stream >> temp; 
			input_stream >> inputParam[temp].asDouble;		// courant_factor
			input_stream >> temp; 
			input_stream >> inputParam[temp].asInt;			// timesteps
			input_stream >> temp;
			input_stream >> inputParam[temp].asInt;			// output_every
			input_stream >> temp;
			input_stream >> inputParam[temp].asInt;			// gridpoints
			input_stream >> temp;
			input_stream >> inputParam[temp].asDouble;		// time_step_size
			input_stream >> temp;
			input_stream >> inputParam[temp].asDouble;		// init_cond_scale
			input_stream >> temp;
			input_stream >> inputParam[temp].asDouble;		// init_cond_width
			input_stream >> temp;
			input_stream >> inputParam[temp].asDouble;		// init_cond_mean
			input_stream >> temp;
			input_stream >> inputParam[temp].asDouble;		// init_cond_sigma

			input_stream.close();
		}
		else
		{
			std::cout << "Cannot open parameter file!" << std::endl;
			throw 0;
		}
	}
	/****************** END: read_parameter *********************/
	
	simple_header::simple_header (int i_n_data_headers, std::string *i_data_headers, std::string i_comment_string, int i_logger) : header (i_logger) {
		int i;

		TRACE (logger, "Constructing...");

		header_string = i_comment_string;
			
		for (i = 0; i < i_n_data_headers; ++i)
		{
			header_string.append (" ");
			header_string.append (i_data_headers [i]);
		}
		
		header_string.append ("\n");

		TRACE (logger, "Constructed.");
	}
	
	output::output (header *i_header_ptr, int i_n, int i_logger) {
		logger = i_logger;
		
		TRACE (logger, "Instantiating...");
		
		n = i_n;
		n_data_ptrs = 0;
		header_ptr = i_header_ptr;

		TRACE (logger, "Instantiated.");
	}
	
	void output::append (double *data_ptr) {
		TRACE (logger, "Appending double to output...");
		
		++n_data_ptrs;
		double_ptrs.push_back (data_ptr);
		int_ptrs.push_back (NULL);
		
		TRACE (logger, "Appended.");
	}
	
	void output::append (int *data_ptr) {
		TRACE (logger, "Appending int to output...");
		
		++n_data_ptrs;
		double_ptrs.push_back (NULL);
		int_ptrs.push_back (data_ptr);
		
		TRACE (logger, "Appended.");
	}
	
	void output::std_to_file (std::string file_name) {
		int i, j;
		std::ofstream file_stream; // A file stream object to be used when writing to file 
		
		TRACE (logger, "Beginning output...");	
				
		INFO (logger, "Outputting to file " << file_name << "...");
		
		file_stream.open (file_name.c_str ());
		
		if (! file_stream.is_open ()) {
			exceptions::file_exception failure;
			ERROR (logger, "Failed to open file " << file_name);
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
		
		TRACE (logger, "Output to file.");
	}
	
	std::string incremental_output::generate_file_name () {

		TRACE (logger, "Generating new file_name...");
		
		if (output_count / output_every >= std::pow (10., int_width)) {
			ERROR (logger, "Exceeded maximum number of outputs")
		}
		std::ostringstream new_file_name;
		
		new_file_name << file_base << std::setfill ('0') << std::setw (int_width) << output_count / output_every << file_extension;

		TRACE (logger, "New file name is " << new_file_name.str () << ".");
		
		return std::string (new_file_name.str ());
	}
} /* io */
