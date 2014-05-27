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
	// netCDF::NcError error_instance (netCDF::NcError::verbose_nonfatal);
	std::map <std::string, virtual_dump> virtual_dumps;
	
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
	
	namespace one_d
	{
		// void ascii::open_file (std::string file_name, int file_type) {
		// 	if (file_type == format::read_file) {
		// 		file_stream.open (file_name, std::ostream::in);
		// 	} else if (file_type == format::replace_file) {
		// 		file_stream.open (file_name, std::ostream::out | std::ostream::trunc);
		// 	} else {
		// 		throw 0;
		// 	}
		// 	
		// 	double_data.resize (0);
		// 	int_data.resize (0);
		// 	float_data.resize (0);
		// 	types.resize (0);
		// 	names.resize (0);
		// 	
		// 	if (! file_stream.is_open ()) {
		// 		ERROR ("Failed to open file " << file_name);
		// 		throw exceptions::file_exception ();
		// 	}
		// }
		// 
		// void ascii::close_file () {
		// 	TRACE ("Closing file.");
		// 	file_stream << comment;
		// 	for (int j = 0; j < (int) names.size (); ++j) {
		// 		file_stream << ' ' << names [j];
		// 	}
		// 	file_stream << '\n';
		// 	
		// 	for (int i = 0; i < n; ++i)
		// 	{
		// 		for (int j = 0; j < (int) names.size (); ++j)
		// 		{
		// 			if (*types [j] == typeid (int)) {
		// 				file_stream << int_data [j] [i];
		// 			} else if (*types [j] == typeid (double)) {
		// 				file_stream << double_data [j] [i];
		// 			} else if (*types [j] == typeid (float)) {
		// 				file_stream << float_data [j] [i];
		// 			}
		// 			file_stream << ' ';
		// 		}
		// 		file_stream << '\n';
		// 	}
		// 	
		// 	file_stream.close ();
		// }
		
		// void netcdf::open_file (std::string filename, int file_type) {
		// 	if (file_type == format::read_file) {
		// 		datafile = new netCDF::NcFile (filename.c_str (), netCDF::NcFile::read);
		// 	} else if (file_type & format::replace_file) {
		// 		datafile = new netCDF::NcFile (filename.c_str (), netCDF::NcFile::replace);
		// 	} else {
		// 		throw 0;
		// 	}
		// 				
		// 	// if (!(((netCDF::NcFile *) datafile)->is_valid ())) {
		// 	// 	FATAL ("Unable to output.");
		// 	// 	throw 0;
		// 	// }
		// 	
		// 	zDim = ((netCDF::NcFile *) datafile)->addDim ("z", n);
		// }
		// 
		// void netcdf::close_file () {
		// 	delete ((netCDF::NcFile *) datafile);
		// }
		// 
		// void netcdf::write (std::string name, int *data) {
		// 	netCDF::NcVar ncdata = ((netCDF::NcFile *) datafile)->addVar (name.c_str (), netCDF::ncInt, zDim);
		// 	ncdata.set_cur (n_offset);
		// 	ncdata.putVar ((int *) (data), n);
		// }
		// void netcdf::write (std::string name, double *data) {
		// 	netCDF::NcVar ncdata = ((netCDF::NcFile *) datafile)->addVar (name.c_str (), netCDF::ncDouble, zDim);
		// 	ncdata.set_cur (n_offset);
		// 	ncdata.putVar ((double *) (data), n);
		// }
		// void netcdf::write (std::string name, float *data) {
		// 	netCDF::NcVar ncdata = ((netCDF::NcFile *) datafile)->addVar (name.c_str (), netCDF::ncFloat, zDim);
		// 	ncdata.set_cur (n_offset);
		// 	ncdata.putVar ((float *) (data), n);
		// }
		// 
		// void netcdf::write_scalar (std::string name, int *data) {
		// 	netCDF::NcVar ncdata = ((netCDF::NcFile *) datafile)->addVar (name.c_str (), netCDF::ncInt);
		// 	ncdata.putVar ((int *) (data));
		// }
		// void netcdf::write_scalar (std::string name, double *data) {
		// 	netCDF::NcVar ncdata = ((netCDF::NcFile *) datafile)->addVar (name.c_str (), netCDF::ncDouble);
		// 	ncdata.putVar ((double *) (data));
		// }
		// void netcdf::write_scalar (std::string name, float *data) {
		// 	netCDF::NcVar ncdata = ((netCDF::NcFile *) datafile)->addVar (name.c_str (), netCDF::ncFloat);
		// 	ncdata.putVar ((float *) (data));
		// }
		// 
		// void netcdf::read (std::string name, int *data) {
		// 	netCDF::NcVar ncdata = ((netCDF::NcFile *) datafile)->get_var (name.c_str ());
		// 	ncdata.set_cur (n_offset);
		// 	ncdata.get ((int *) (data), n);
		// }
		// void netcdf::read (std::string name, double *data) {
		// 	netCDF::NcVar ncdata = ((netCDF::NcFile *) datafile)->get_var (name.c_str ());
		// 	ncdata.set_cur (n_offset);
		// 	ncdata.get ((double *) (data), n);
		// }
		// void netcdf::read (std::string name, float *data) {
		// 	netCDF::NcVar ncdata = ((netCDF::NcFile *) datafile)->get_var (name.c_str ());
		// 	ncdata.set_cur (n_offset);
		// 	ncdata.get ((float *) (data), n);
		// }
		// 
		// void netcdf::read_scalar (std::string name, int *data) {
		// 	netCDF::NcVar ncdata = ((netCDF::NcFile *) datafile)->get_var (name.c_str ());
		// 	ncdata.get ((int *) (data));
		// }
		// void netcdf::read_scalar (std::string name, double *data) {
		// 	netCDF::NcVar ncdata = ((netCDF::NcFile *) datafile)->get_var (name.c_str ());
		// 	ncdata.get ((double *) (data));
		// }
		// void netcdf::read_scalar (std::string name, float *data) {
		// 	netCDF::NcVar ncdata = ((netCDF::NcFile *) datafile)->get_var (name.c_str ());
		// 	ncdata.get ((float *) (data));
		// }
	} /* one_d */
	
	namespace two_d
	{
		std::map <std::string, netCDF::NcFile *> netcdf::files;
		std::map <std::string, std::vector <netCDF::NcDim>> netcdf::dims;
		std::map <std::string, std::vector <std::string>> netcdf::failures;
		
		void netcdf::open_file (std::string file_name, int file_type, int n_max, int m_max, int l_max) {
			if (files [file_name]) {
				FATAL ("File not closed.");
				throw 0;
			}
			if (file_type == read_file) {
				files [file_name] = new netCDF::NcFile (file_name.c_str (), netCDF::NcFile::read);
			} else if (file_type == replace_file) {
				files [file_name] = new netCDF::NcFile (file_name.c_str (), netCDF::NcFile::replace);
				dims [file_name].resize (2);
				dims [file_name] [0] = files [file_name]->addDim ("x", n_max);
				dims [file_name] [1] = files [file_name]->addDim ("z", m_max);
			} else {
				throw 0;
			}
			
			if (files [file_name]->isNull ()) {
				FATAL ("Unable to load file.");
				throw 0;
			}
			
			failures [file_name].resize (0);
		}
		
		void netcdf::close_file (std::string file_name) {
			
			delete files [file_name];
			
			if (failures [file_name].size () > 0) {
				throw exceptions::io::bad_variables ((int) failures [file_name].size (), &failures [file_name] [0]);
			}
			
			/*
				TODO Move exception to input class
			*/
		}
	} /* two_d */
} /* io */
