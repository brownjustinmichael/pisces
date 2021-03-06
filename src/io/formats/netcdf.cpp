/*!**********************************************************************
 * \file netcdf.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "netcdf.hpp"

namespace formats
{
	std::map <std::string, netCDF::NcFile *> netcdf::files;
	std::map <std::string, std::vector <netCDF::NcDim>> netcdf::dims;
	std::map <std::string, std::vector <std::string>> netcdf::failures;
	std::map <std::string, int> netcdf::records;
	std::map <std::string, bool> netcdf::first;
	bool netcdf::uses_files = true;

	void netcdf::open_file (const data_grid &grid, std::string file_name, int file_type) {
		TRACE ("Opening NetCDF file");
		// Update records and first
		if (!first [file_name]) {
			records [file_name] = 0;
			first [file_name] = true;
		} else {
			records [file_name] = -1;
		}
		
		if (file_type == read_file) {
			if (files [file_name]) {
				FATAL ("File already open");
				throw 0;
			}
			// Load information from file
			files [file_name] = new netCDF::NcFile (file_name.c_str (), netCDF::NcFile::read);
			if (files [file_name]->isNull () or files [file_name]->getDim ("time").isNull ()) {
				records [file_name] = 0;
			} else {
				records [file_name] = files [file_name]->getDim ("time").getSize () - 1;
			}
		} else if (file_type == replace_file || records [file_name] == 0) {
			if (files [file_name]) {
				FATAL ("File already open");
				throw 0;
			}
			// Set up dimensions
			files [file_name] = new netCDF::NcFile (file_name.c_str (), netCDF::NcFile::replace);
			dims [file_name].resize (3);
			dims [file_name] [0] = files [file_name]->addDim ("time");
			dims [file_name] [1] = files [file_name]->addDim ("x", grid.get_max (0));
			dims [file_name] [2] = files [file_name]->addDim ("z", grid.get_max (1));
			records [file_name] = 0.;
		} else if (file_type == append_file) {
			// Set up dimensions
			files [file_name] = new netCDF::NcFile (file_name.c_str (), netCDF::NcFile::write);
			dims [file_name] [0] = files [file_name]->getDim ("time");
			dims [file_name] [1] = files [file_name]->getDim ("x");
			dims [file_name] [2] = files [file_name]->getDim ("z");
			if (files [file_name]->isNull () or files [file_name]->getDim ("time").isNull ()) {
				records [file_name] = 0;
			} else {
				records [file_name] = files [file_name]->getDim ("time").getSize ();
			}
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
		records [file_name] += 1;
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
} /* formats */
