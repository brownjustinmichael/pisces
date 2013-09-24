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
			input_stream >> inputParam[temp].asInt;			// n_elements
			input_stream >> temp;
			input_stream >> inputParam[temp].asInt;			// n_iterations
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
	
	simple_header::simple_header (int i_n_data_headers, std::string *i_data_headers, std::string i_comment_string) {
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
	
	template <class datatype>
	output <datatype>::output (header *i_header_ptr, int i_n) {
		TRACE ("Instantiating...");
		
		n = i_n;
		n_data_ptrs = 0;
		header_ptr.reset (i_header_ptr);

		TRACE ("Instantiated.");
	}
	
	template <class datatype>
	void output <datatype>::append (datatype *data_ptr) {
		TRACE ("Appending datatype to output...");
		
		++n_data_ptrs;
		datatype_ptrs.push_back (data_ptr);
		int_ptrs.push_back (NULL);
		
		TRACE ("Appended.");
	}
	
	template <class datatype>
	void output <datatype>::append (int *data_ptr) {
		TRACE ("Appending int to output...");
		
		++n_data_ptrs;
		datatype_ptrs.push_back (NULL);
		int_ptrs.push_back (data_ptr);
		
		TRACE ("Appended.");
	}
	
	template <class datatype>
	void output <datatype>::std_to_file (std::string file_name) {
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
				if (datatype_ptrs [j] == NULL) {
					file_stream << int_ptrs [j] [i];
				} else {
					file_stream << datatype_ptrs [j] [i];
				}
				file_stream << ' ';
			}
			file_stream << '\n';
		}
				
		file_stream.close ();
		
		TRACE ("Output to file.");
	}
	
	template <class datatype>
	std::string incremental_output <datatype>::generate_file_name () {

		TRACE ("Generating new file_name...");
		
		if (output_count / output_every >= std::pow (10., int_width)) {
			ERROR ("Exceeded maximum number of outputs")
		}
		std::ostringstream new_file_name;
		
		new_file_name << file_base << std::setfill ('0') << std::setw (int_width) << output_count / output_every << file_extension;

		TRACE ("New file name is " << new_file_name.str () << ".");
		
		return std::string (new_file_name.str ());
	}
	
	template class output <double>;
	template class output <float>;
	
	template class simple_output <double>;
	template class simple_output <float>;
	
	template class incremental_output <double>;
	template class incremental_output <float>;
} /* io */
