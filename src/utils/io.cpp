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
	
	namespace one_d
	{
		void ascii::open_file (std::string file_name, int file_type) {
			if (file_type == format::read_file) {
				file_stream.open (file_name, std::ostream::in);
			} else if (file_type == format::replace_file) {
				file_stream.open (file_name, std::ostream::out | std::ostream::trunc);
			} else {
				throw 0;
			}
			
			double_data.resize (0);
			int_data.resize (0);
			float_data.resize (0);
			types.resize (0);
			names.resize (0);
			
			if (! file_stream.is_open ()) {
				ERROR ("Failed to open file " << file_name);
				throw exceptions::file_exception ();
			}
		}
		
		void ascii::close_file () {
			TRACE ("Closing file.");
			file_stream << comment;
			for (int j = 0; j < (int) names.size (); ++j) {
				file_stream << ' ' << names [j];
			}
			file_stream << '\n';
			
			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j < (int) names.size (); ++j)
				{
					if (*types [j] == typeid (int)) {
						file_stream << int_data [j] [i];
					} else if (*types [j] == typeid (double)) {
						file_stream << double_data [j] [i];
					} else if (*types [j] == typeid (float)) {
						file_stream << float_data [j] [i];
					}
					file_stream << ' ';
				}
				file_stream << '\n';
			}
			
			file_stream.close ();
		}
		
		void ascii::write (std::string name, double *data) {
			TRACE ("Writing...");
			names.push_back (name);
			double_data.resize (double_data.size () + 1);
			float_data.resize (float_data.size () + 1);
			int_data.resize (int_data.size () + 1);
			double_data [(int) double_data.size () - 1].resize (n);
			double_data [(int) double_data.size () - 1].assign (data, data + n);
			types.push_back (&typeid (double));
		}
		void ascii::write (std::string name, float *data) {
			TRACE ("Writing...");
			names.push_back (name);
			double_data.resize (double_data.size () + 1);
			float_data.resize (float_data.size () + 1);
			int_data.resize (int_data.size () + 1);
			float_data [(int) float_data.size () - 1].resize (n);
			float_data [(int) float_data.size () - 1].assign (data, data + n);
			types.push_back (&typeid (float));
		}
		void ascii::write (std::string name, int *data) {
			TRACE ("Writing...");
			names.push_back (name);
			double_data.resize (double_data.size () + 1);
			float_data.resize (float_data.size () + 1);
			int_data.resize (int_data.size () + 1);
			int_data [(int) int_data.size () - 1].resize (n);
			int_data [(int) int_data.size () - 1].assign (data, data + n);
			types.push_back (&typeid (int));
		}
		
		void ascii::write_scalar (std::string name, double *data) {
		}
		void ascii::write_scalar (std::string name, float *data) {
		}
		void ascii::write_scalar (std::string name, int *data) {
		}
		
		void ascii::read (std::string name, double *data) {
			throw 0;
		}
		void ascii::read (std::string name, float *data) {
			throw 0;
		}
		void ascii::read (std::string name, int *data) {
			throw 0;
		}
		
		void ascii::read_scalar (std::string name, double *data) {
			throw 0;
		}
		void ascii::read_scalar (std::string name, float *data) {
			throw 0;
		}
		void ascii::read_scalar (std::string name, int *data) {
			throw 0;
		}
		
		void netcdf::open_file (std::string filename, int file_type) {
			if (file_type == format::read_file) {
				datafile = new NcFile (filename.c_str (), NcFile::ReadOnly);
			} else if (file_type & format::replace_file) {
				datafile = new NcFile (filename.c_str (), NcFile::Replace);
			} else {
				throw 0;
			}
						
			if (!(((NcFile *) datafile)->is_valid ())) {
				FATAL ("Unable to output.");
				throw 0;
			}
			
			zDim = ((NcFile *) datafile)->add_dim ("z", n);
		}
		
		void netcdf::close_file () {
			delete ((NcFile *) datafile);
		}
		
		void netcdf::write (std::string name, int *data) {
			NcVar *ncdata = ((NcFile *) datafile)->add_var (name.c_str (), ncInt, (NcDim *) zDim);
			ncdata->put ((int *) (data), n);
		}
		void netcdf::write (std::string name, double *data) {
			NcVar *ncdata = ((NcFile *) datafile)->add_var (name.c_str (), ncDouble, (NcDim *) zDim);
			ncdata->put ((double *) (data), n);
		}
		void netcdf::write (std::string name, float *data) {
			NcVar *ncdata = ((NcFile *) datafile)->add_var (name.c_str (), ncFloat, (NcDim *) zDim);
			ncdata->put ((float *) (data), n);
		}
		
		void netcdf::write_scalar (std::string name, int *data) {
			NcVar *ncdata = ((NcFile *) datafile)->add_var (name.c_str (), ncInt);
			ncdata->put ((int *) (data));
		}
		void netcdf::write_scalar (std::string name, double *data) {
			NcVar *ncdata = ((NcFile *) datafile)->add_var (name.c_str (), ncDouble);
			ncdata->put ((double *) (data));
		}
		void netcdf::write_scalar (std::string name, float *data) {
			NcVar *ncdata = ((NcFile *) datafile)->add_var (name.c_str (), ncFloat);
			ncdata->put ((float *) (data));
		}
		
		void netcdf::read (std::string name, int *data) {
			NcVar *ncdata = ((NcFile *) datafile)->get_var (name.c_str ());
			ncdata->get ((int *) (data), n);
		}
		void netcdf::read (std::string name, double *data) {
			NcVar *ncdata = ((NcFile *) datafile)->get_var (name.c_str ());
			ncdata->get ((double *) (data), n);
		}
		void netcdf::read (std::string name, float *data) {
			NcVar *ncdata = ((NcFile *) datafile)->get_var (name.c_str ());
			ncdata->get ((float *) (data), n);
		}
		
		void netcdf::read_scalar (std::string name, int *data) {
			NcVar *ncdata = ((NcFile *) datafile)->get_var (name.c_str ());
			ncdata->get ((int *) (data));
		}
		void netcdf::read_scalar (std::string name, double *data) {
			NcVar *ncdata = ((NcFile *) datafile)->get_var (name.c_str ());
			ncdata->get ((double *) (data));
		}
		void netcdf::read_scalar (std::string name, float *data) {
			NcVar *ncdata = ((NcFile *) datafile)->get_var (name.c_str ());
			ncdata->get ((float *) (data));
		}
	} /* one_d */
	
	namespace two_d
	{
		void netcdf::open_file (std::string filename, int file_type) {
			if (file_type == format::read_file) {
				datafile = new NcFile (filename.c_str (), NcFile::ReadOnly);
			} else if (file_type & format::replace_file) {
				datafile = new NcFile (filename.c_str (), NcFile::Replace);
			} else {
				throw 0;
			}
			
			if (!(((NcFile *) datafile)->is_valid ())) {
				FATAL ("Unable to load file.");
				throw 0;
			}
			
			xDim = ((NcFile *) datafile)->add_dim ("x", n);
			zDim = ((NcFile *) datafile)->add_dim ("z", m);
			
			failures.resize (0);
		}
		
		void netcdf::close_file () {
			delete ((NcFile *) datafile);
			
			if (failures.size () > 0) {
				throw exceptions::io::bad_variables ((int) failures.size (), &failures [0]);
			}
			
			/*
				TODO Move exception to input class
			*/
		}
		
		void netcdf::write (std::string name, int *data) {
			NcVar *ncdata = ((NcFile *) datafile)->add_var (name.c_str (), ncInt, (NcDim *) xDim, (NcDim *) zDim);
			ncdata->put ((int *) (data), n, m);
		}
		void netcdf::write (std::string name, double *data) {
			NcVar *ncdata = ((NcFile *) datafile)->add_var (name.c_str (), ncDouble, (NcDim *) xDim, (NcDim *) zDim);
			ncdata->put ((double *) (data), n, m);
		}
		void netcdf::write (std::string name, float *data) {
			NcVar *ncdata = ((NcFile *) datafile)->add_var (name.c_str (), ncFloat, (NcDim *) xDim, (NcDim *) zDim);
			ncdata->put ((float *) (data), n, m);
		}
		
		void netcdf::write_scalar (std::string name, int *data) {
			NcVar *ncdata = ((NcFile *) datafile)->add_var (name.c_str (), ncInt);
			ncdata->put ((int *) (data));
		}
		void netcdf::write_scalar (std::string name, double *data) {
			NcVar *ncdata = ((NcFile *) datafile)->add_var (name.c_str (), ncDouble);
			ncdata->put ((double *) (data));
		}
		void netcdf::write_scalar (std::string name, float *data) {
			NcVar *ncdata = ((NcFile *) datafile)->add_var (name.c_str (), ncFloat);
			ncdata->put ((float *) (data));
		}
		
		void netcdf::read (std::string name, int *data) {
			NcVar *ncdata = ((NcFile *) datafile)->get_var (name.c_str ());
			if (ncdata && ncdata->is_valid ()) {
				ncdata->get ((int *) (data), n, m);
			} else {
				failures.push_back (name);
				WARN ("Variable " << name << " not found in file");
			}
		}
		void netcdf::read (std::string name, double *data) {
			NcVar *ncdata = ((NcFile *) datafile)->get_var (name.c_str ());
			if (ncdata && ncdata->is_valid ()) {
				ncdata->get ((double *) (data), n, m);
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j < m; ++j) {
						if (std::isnan (data [i * m + j])) {
							WARN ("Reading in Nan.");
							throw 0;
						}
					}
				}
			} else {
				failures.push_back (name);
				WARN ("Variable " << name << " not found in file");
			}
		}
		void netcdf::read (std::string name, float *data) {
			NcVar *ncdata = ((NcFile *) datafile)->get_var (name.c_str ());
			if (ncdata && ncdata->is_valid ()) {
				ncdata->get ((float *) (data), n, m);
			} else {
				failures.push_back (name);
				WARN ("Variable " << name << " not found in file");
			}
		}
		
		void netcdf::read_scalar (std::string name, int *data) {
			NcVar *ncdata = ((NcFile *) datafile)->get_var (name.c_str ());
			if (ncdata && ncdata->is_valid ()) {
				ncdata->get ((int *) (data));
			} else {
				failures.push_back (name);
				WARN ("Variable " << name << " not found in file");
			}
		}
		void netcdf::read_scalar (std::string name, double *data) {
			NcVar *ncdata = ((NcFile *) datafile)->get_var (name.c_str ());
			if (ncdata && ncdata->is_valid ()) {
				ncdata->get ((double *) (data));
			} else {
				failures.push_back (name);
				WARN ("Variable " << name << " not found in file");
			}
		}
		void netcdf::read_scalar (std::string name, float *data) {
			NcVar *ncdata = ((NcFile *) datafile)->get_var (name.c_str ());
			if (ncdata && ncdata->is_valid ()) {
				ncdata->get ((float *) (data));
			} else {
				failures.push_back (name);
				WARN ("Variable " << name << " not found in file");
			}
		}
	} /* two_d */
} /* io */
