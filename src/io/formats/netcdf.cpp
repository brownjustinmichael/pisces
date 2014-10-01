/*!**********************************************************************
 * \file netcdf.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "netcdf.hpp"

namespace io
{
	namespace formats
	{
		namespace two_d
		{
			std::map <std::string, netCDF::NcFile *> netcdf::files;
			std::map <std::string, std::vector <netCDF::NcDim>> netcdf::dims;
			std::map <std::string, std::vector <std::string>> netcdf::failures;
			std::map <std::string, int> netcdf::records;
		
			void netcdf::open_file (std::string file_name, int file_type, int n_max, int m_max, int l_max) {
				if (file_type == read_file) {
					if (files [file_name]) {
						FATAL ("File already open");
						throw 0;
					}
					files [file_name] = new netCDF::NcFile (file_name.c_str (), netCDF::NcFile::read);
					if (files [file_name]->isNull () or files [file_name]->getDim ("record").isNull ()) {
						records [file_name] = 0;
					} else {
						records [file_name] = files [file_name]->getDim ("record").getSize () - 1;
					}
				} else if (file_type == replace_file) {
					if (files [file_name]) {
						FATAL ("File already open");
						throw 0;
					}
					files [file_name] = new netCDF::NcFile (file_name.c_str (), netCDF::NcFile::replace);
					dims [file_name].resize (3);
					dims [file_name] [0] = files [file_name]->addDim ("record");
					dims [file_name] [1] = files [file_name]->addDim ("x", n_max);
					dims [file_name] [2] = files [file_name]->addDim ("z", m_max);
					records [file_name] = 0;
				} else if (file_type == append_file) {
					dims [file_name].resize (3);
					if (!files [file_name]) {
						records [file_name] = 0;
						files [file_name] = new netCDF::NcFile (file_name.c_str (), netCDF::NcFile::replace);
						dims [file_name] [0] = files [file_name]->addDim ("record");
						dims [file_name] [1] = files [file_name]->addDim ("x", n_max);
						dims [file_name] [2] = files [file_name]->addDim ("z", m_max);
					}
					records [file_name] = dims [file_name] [0].getSize ();
				} else {
					FATAL ("Unrecognized file type.");
					throw 0;
				}
			
				if (files [file_name]->isNull ()) {
					FATAL ("Unable to load file.");
					throw 0;
				}
			
				failures [file_name].resize (0);
			}
		
			void netcdf::close_file (std::string file_name, int file_type) {
				if (files [file_name]) {
					delete files [file_name];
					files [file_name] = NULL;
				}
			
			
				if (failures [file_name].size () > 0) {
					// throw exceptions::io::bad_variables ();
				}
			
				/*
					TODO Move exception to input class
				*/
			}
		} /* two_d */
	} /* formats */
} /* io */