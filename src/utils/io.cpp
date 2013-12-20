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
	template <class datatype>
	parameters <datatype>::parameters (std::string i_filename) {
		std::string word;
		std::ifstream input_stream (i_filename);
		std::map <std::string, std::string> read_map;
		
		while (!(input_stream.eof ())) {
			input_stream >> word;
			input_stream >> read_map [word];
		}
		
		diffusion_coeff = atof (read_map ["diffusion_coeff"].c_str ());
		nonlinear_diffusion_coeff = atof (read_map ["nonlinear_diffusion_coeff"].c_str ());
		advection_coeff = atof (read_map ["advection_coeff"].c_str ());
		courant_factor = atof (read_map ["courant_factor"].c_str ());
		timesteps = atoi (read_map ["timesteps"].c_str ());
		output_every = atoi (read_map ["output_every"].c_str ());
		gridpoints = atoi (read_map ["gridpoints"].c_str ());
		n_iterations = atoi (read_map ["n_iterations"].c_str ());
		max_timestep = atof (read_map ["max_timestep"].c_str ());
		scale = atof (read_map ["scale"].c_str ());
		mean = atof (read_map ["mean"].c_str ());
		width = atof (read_map ["width"].c_str ());
		sigma = atof (read_map ["sigma"].c_str ());
		implicit_alpha = atof (read_map ["implicit_alpha"].c_str ());
		n = atoi (read_map ["n"].c_str ());
		nrhs = atoi (read_map ["nrhs"].c_str ());
		nmp = atoi (read_map ["nmp"].c_str ());
		nb = atoi (read_map ["nb"].c_str ());
		output = read_map ["output"].c_str ();
	}
	
	template class parameters <double>;
	template class parameters <float>;
	
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
