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
		
		while (!(input_stream.eof ())) {
			input_stream >> word;
			if (word == "diffusion_coeff") {
				input_stream >> diffusion_coeff;
			} else if (word == "advection_coeff") {
				input_stream >> advection_coeff;
			} else if (word == "courant_factor") {
				input_stream >> courant_factor;
			} else if (word == "timesteps") {
				input_stream >> timesteps;
			} else if (word == "output_every") {
				input_stream >> output_every;
			} else if (word == "gridpoints") {
				input_stream >> gridpoints;
			} else if (word == "n_iterations") {
				input_stream >> n_iterations;
			} else if (word == "max_timestep") {
				input_stream >> max_timestep;
			} else if (word == "scale") {
				input_stream >> scale;
			} else if (word == "width") {
				input_stream >> width;
			} else if (word == "mean") {
				input_stream >> mean;
			} else if (word == "sigma") {
				input_stream >> sigma;
			} else if (word == "n") {
				input_stream >> n;
			} else if (word == "nrhs") {
				input_stream >> nrhs;
			} else if (word == "nmp") {
				input_stream >> nmp;
			} else if (word == "nb") {
				input_stream >> nb;
			}
		}
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
