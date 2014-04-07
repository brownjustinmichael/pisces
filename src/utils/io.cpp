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
	NcError error_instance (NcError::verbose_nonfatal);
	
	YAML::Node parameters::operator [] (std::string key) {
		std::istringstream ss (key);
		std::string token;
		std::getline (ss, token, '.');
		std::vector <YAML::Node> nodes;
		nodes.push_back (YAML::Node::operator [] (token));
		while (std::getline (ss, token, '.')) {
			nodes.push_back (nodes [(int) nodes.size () - 1] [token]);
		}
		return nodes [(int) nodes.size () - 1];
	}
	
	template <>
	std::map <std::string, std::vector <double>>::iterator virtual_dump::begin <double> () {
		return double_map.begin ();
	}
	
	template <>
	std::map <std::string, std::vector <float>>::iterator virtual_dump::begin <float> () {
		return float_map.begin ();
	}
	
	template <>
	std::map <std::string, std::vector <int>>::iterator virtual_dump::begin <int> () {
		return int_map.begin ();
	}
	
	template <>
	std::map <std::string, std::vector <double>>::iterator virtual_dump::end <double> () {
		return double_map.end ();
	}
	
	template <>
	std::map <std::string, std::vector <float>>::iterator virtual_dump::end <float> () {
		return float_map.end ();
	}
	
	template <>
	std::map <std::string, std::vector <int>>::iterator virtual_dump::end <int> () {
		return int_map.end ();
	}
	
	template <>
	double &virtual_dump::index <double> (std::string name, int i, int j) {
		TRACE ("Indexing " << name << " " << &double_map [name] [i]);
		TRACE (" " << dims [name] [0]);
		return double_map [name] [i * dims [name] [1] + j];
		TRACE ("Done.");
	}
	
	template <>
	float &virtual_dump::index <float> (std::string name, int i, int j) {
		TRACE ("Indexing " << name);
		return float_map [name] [i * dims [name] [1] + j];
	}
	
	template <>
	int &virtual_dump::index <int> (std::string name, int i, int j) {
		TRACE ("Indexing " << name);
		return int_map [name] [i * dims [name] [1] + j];
	}
	
	namespace one_d
	{
		void ascii::to_file (const char *file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs) {
			int i, j;
			std::ofstream file_stream; // A file stream object to be used when writing to file 
		
			TRACE ("Beginning output...");	
				
			INFO ("Outputting to file " << file_name << "...");
		
			file_stream.open (file_name);
		
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
		
		void ascii::from_file (const char *file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs) {
			// int i, j;
			std::ofstream file_stream; // A file stream object to be used when writing to file 
		
			TRACE ("Beginning input...");	
				
			INFO ("Inputting from file " << file_name << "...");
		
			file_stream.open (file_name);
		
			if (! file_stream.is_open ()) {
				exceptions::file_exception failure;
				ERROR ("Failed to open file " << file_name);
				throw failure;
			}
							
			// for (i = 0; i < n; ++i)
			// {
			// 	for (j = 0; j < n_data_ptrs; ++j)
			// 	{
			// 		if (*types [j] == typeid (int)) {
			// 			file_stream << ((int *) (data_ptrs [j])) [i];
			// 		} else if (*types [j] == typeid (double)) {
			// 			file_stream << ((double *) (data_ptrs [j])) [i];
			// 		} else if (*types [j] == typeid (float)) {
			// 			file_stream << ((float *) (data_ptrs [j])) [i];
			// 		}
			// 		file_stream << ' ';
			// 	}
			// 	file_stream << '\n';
			// }
			
			/*
				TODO Implement this
			*/
				
			file_stream.close ();
		
			TRACE ("Output to file.");
		}
		
		void netcdf::to_file (const char *file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs) {
			NcFile datafile (file_name, NcFile::Replace);
	
			INFO ("Outputting to file " << file_name << "...");
	
			if (!datafile.is_valid ()) {
				FATAL ("Unable to output.");
				throw 0;
			}
	
			NcDim* zDim = datafile.add_dim ("z", n);
		
			for (int j = 0; j < n_data_ptrs; ++j) {
				if (*types [j] == typeid (int)) {
					NcVar *data = datafile.add_var (names [j].c_str (), ncInt, zDim);
					data->put ((int *) (data_ptrs [j]), n);
				} else if (*types [j] == typeid (double)) {
					NcVar *data = datafile.add_var (names [j].c_str (), ncDouble, zDim);
					data->put ((double *) (data_ptrs [j]), n);
				} else if (*types [j] == typeid (float)) {
					NcVar *data = datafile.add_var (names [j].c_str (), ncFloat, zDim);
					data->put ((float *) (data_ptrs [j]), n);
				}
			}
			
			for (int j = 0; j < n_scalar_ptrs; ++j) {
				if (*scalar_types [j] == typeid (int)) {
					NcVar *data = datafile.add_var (scalar_names [j].c_str (), ncInt);
					data->put ((int *) (scalar_ptrs [j]));
				} else if (*scalar_types [j] == typeid (double)) {
					NcVar *data = datafile.add_var (scalar_names [j].c_str (), ncDouble);
					data->put ((double *) (scalar_ptrs [j]));
				} else if (*scalar_types [j] == typeid (float)) {
					NcVar *data = datafile.add_var (names [j].c_str (), ncFloat);
					data->put ((float *) (scalar_ptrs [j]));
				}
			}
		}
		
		void netcdf::from_file (const char *file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs) {
			NcFile datafile (file_name, NcFile::ReadOnly);
	
			INFO ("Inputting from file " << file_name << "...");
	
			if (!datafile.is_valid ()) {
				FATAL ("Could not open " << file_name);
				throw 0;
			}
	
			for (int j = 0; j < n_data_ptrs; ++j) {
				NcVar *data = datafile.get_var (names [j].c_str ());
				if (*types [j] == typeid (int)) {
					data->get ((int *) (data_ptrs [j]), n);
				} else if (*types [j] == typeid (double)) {
					data->get ((double *) (data_ptrs [j]), n);
				} else if (*types [j] == typeid (float)) {
					data->get ((float *) (data_ptrs [j]), n);
				}
			}
			
			for (int j = 0; j < n_scalar_ptrs; ++j) {
				NcVar *data = datafile.get_var (scalar_names [j].c_str ());
				if (*scalar_types [j] == typeid (int)) {
					data->get ((int *) (scalar_ptrs [j]));
				} else if (*scalar_types [j] == typeid (double)) {
					data->get ((double *) (scalar_ptrs [j]));
				} else if (*scalar_types [j] == typeid (float)) {
					data->get ((float *) (scalar_ptrs [j]));
				}
			}
		}
	} /* one_d */
	
	namespace two_d
	{
		void netcdf::to_file (const char *file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs) {
			NcFile *datafile = new NcFile (file_name, NcFile::Replace);
				
			INFO ("Outputting to file " << file_name << "...");
			
			if (!datafile->is_valid ()) {
				FATAL ("Unable to output.");
				throw 0;
			}
				
			NcDim* xDim = datafile->add_dim ("x", n);
			NcDim* zDim = datafile->add_dim ("z", m);
					
			for (int j = 0; j < n_data_ptrs; ++j) {
				if (*types [j] == typeid (int)) {
					NcVar *data = datafile->add_var (names [j].c_str (), ncInt, xDim, zDim);
					data->put ((int *) (data_ptrs [j]), n, m);
				} else if (*types [j] == typeid (double)) {
					NcVar *data = datafile->add_var (names [j].c_str (), ncDouble, xDim, zDim);
					data->put ((double *) (data_ptrs [j]), n, m);
				} else if (*types [j] == typeid (float)) {
					NcVar *data = datafile->add_var (names [j].c_str (), ncFloat, xDim, zDim);
					data->put ((float *) (data_ptrs [j]), n, m);
				}
			}
			
			for (int j = 0; j < n_scalar_ptrs; ++j) {
				if (*scalar_types [j] == typeid (int)) {
					NcVar *data = datafile->add_var (scalar_names [j].c_str (), ncInt);
					data->put ((int *) (scalar_ptrs [j]));
				} else if (*scalar_types [j] == typeid (double)) {
					NcVar *data = datafile->add_var (scalar_names [j].c_str (), ncDouble);
					data->put ((double *) (scalar_ptrs [j]));
				} else if (*scalar_types [j] == typeid (float)) {
					NcVar *data = datafile->add_var (names [j].c_str (), ncFloat);
					data->put ((float *) (scalar_ptrs [j]));
				}
			}
			
			delete (datafile);
		}
		
		void netcdf::from_file (const char *file_name, int n_data_ptrs, std::string *names, const std::type_info **types, void **data_ptrs, int n_scalar_ptrs, std::string *scalar_names, const std::type_info **scalar_types, void **scalar_ptrs) {
			failures.resize (0);
			NcFile *datafile = new NcFile (file_name, NcFile::ReadOnly);
				
			INFO ("Inputting from file " << file_name << "...");

			if (!datafile->is_valid ()) {
				FATAL ("Could not open " << file_name);
				throw 0;
			}
				
			for (int j = 0; j < n_data_ptrs; ++j) {
				NcVar *data = datafile->get_var (names [j].c_str ());
				if (data && data->is_valid ()) {
					if (*types [j] == typeid (int)) {
						data->get ((int *) (data_ptrs [j]), n, m);
					} else if (*types [j] == typeid (double)) {
						data->get ((double *) (data_ptrs [j]), n, m);
					} else if (*types [j] == typeid (float)) {
						data->get ((float *) (data_ptrs [j]), n, m);
					}
				} else {
					failures.push_back (names [j]);
					WARN ("Variable " << names [j] << " not found in file " << file_name);
				}
			}
			
			for (int j = 0; j < n_scalar_ptrs; ++j) {
				NcVar *data = datafile->get_var (scalar_names [j].c_str ());
				if (data->is_valid ()) {
					if (*scalar_types [j] == typeid (int)) {
						data->get ((int *) (scalar_ptrs [j]));
					} else if (*scalar_types [j] == typeid (double)) {
						data->get ((double *) (scalar_ptrs [j]));
					} else if (*scalar_types [j] == typeid (float)) {
						data->get ((float *) (scalar_ptrs [j]));
					}
				} else {
					failures.push_back (names [j]);
					WARN ("Variable " << names [j] << " not found in file " << file_name);
				}
			}

			delete (datafile);
			
			if (failures.size () > 0) {
				throw exceptions::io::bad_variables ((int) failures.size (), &failures [0]);
			}
		}
	} /* two_d */
} /* io */
