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
#include <netcdfcpp.h>
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
	
	namespace one_d
	{
		void ascii::to_file (std::string file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs) {
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
							
			for (i = 0; i < n; ++i)
			{
				for (j = 0; j < n_data_ptrs; ++j)
				{
					if (*types [j] == typeid (int)) {
						file_stream << ((int *) (data_ptrs [j])) [i];
					} else if (*types [j] == typeid (double)) {
						file_stream << ((double *) (data_ptrs [j])) [i];
					} else if (*types [j] == typeid (float)) {
						file_stream << ((float *) (data_ptrs [j])) [i];
					}
					file_stream << ' ';
				}
				file_stream << '\n';
			}
				
			file_stream.close ();
		
			TRACE ("Output to file.");
		}
	} /* one_d */
	
	namespace two_d
	{
		void netcdf::to_file (std::string file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs) {
			NcFile datafile (file_name.c_str (), NcFile::Replace);
	
			INFO ("Outputting to file " << file_name << "...");

			if (!datafile.is_valid ()) {
				FATAL ("Unable to output.");
				throw 0;
			}
	
			NcDim* xDim = datafile.add_dim ("x", n);
			NcDim* zDim = datafile.add_dim ("z", m);
		
			for (int j = 0; j < n_data_ptrs; ++j) {
				if (*types [j] == typeid (int)) {
					NcVar *data = datafile.add_var (names [j].c_str (), ncInt, xDim, zDim);
					data->put ((int *) (data_ptrs [j]), n, m);
				} else if (*types [j] == typeid (double)) {
					NcVar *data = datafile.add_var (names [j].c_str (), ncDouble, xDim, zDim);
					data->put ((double *) (data_ptrs [j]), n, m);
				} else if (*types [j] == typeid (float)) {
					NcVar *data = datafile.add_var (names [j].c_str (), ncFloat, xDim, zDim);
					data->put ((float *) (data_ptrs [j]), n, m);
				}
			}
		}
	} /* two_d */
} /* io */
