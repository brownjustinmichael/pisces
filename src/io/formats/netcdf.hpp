/*!**********************************************************************
 * \file netcdf.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef NETCDF_HPP_72C14A02
#define NETCDF_HPP_72C14A02

#include <netcdf>

#include "logger/logger.hpp"

#include "format.hpp"
#include "exceptions.hpp"

namespace io
{
	namespace formats
	{
		/*!**********************************************************************
		 * \brief Returns the netCDF object associated with a given type
		 * 
		 * \param type A pointer to a type_info object associated with the object type in question
		 * 
		 * This function throws an io::bad_type exception if the type is not known
		 * 
		 * \return The netCDF::NcType object associated with the input type
		 ************************************************************************/
		inline netCDF::NcType netcdf_type (const std::type_info* type) {
			if (type == &typeid (double)) {
				return netCDF::ncDouble;
			} else if (type == &typeid (int)) {
				return netCDF::ncInt;
			} else if (type == &typeid (float)) {
				return netCDF::ncFloat;
			} else {
				FATAL ("Unrecognized NetCDF type");
				throw exceptions::bad_type ();
			}
		}
		
		namespace two_d
		{
			class netcdf
			{
			public:
				netcdf () {}
		
				virtual ~netcdf () {}
		
				static std::string extension () {return ".cdf";}
			
				static void open_file (const data_grid &grid, std::string file_name, int file_type);
		
				static void close_file (std::string file_name, int file_type);
			
				template <class datatype>
				static void write (const data_grid &grid, std::string file_name, std::string name, void *data, int record = -1, int flags = all_d) {
					std::vector <netCDF::NcDim> scalar_dims;
					std::vector <size_t> offsets = grid.offsets;
					std::vector <size_t> sizes = grid.ns;
					if (flags < 0) {
						scalar_dims = dims [file_name];
					} else {
						scalar_dims = {dims [file_name] [0]};
						offsets = {offsets [0]};
						sizes = {sizes [0]};
						int i = 0;
						while (flags > 0) {
							if (flags % 2 != 0) scalar_dims.push_back (dims [file_name] [i + 1]);
							if (flags % 2 != 0) offsets.push_back (grid.offsets [i + 1]);
							if (flags % 2 != 0) sizes.push_back (grid.ns [i + 1]);
							flags = flags >> 1;
							i++;
						}
					}
					
					offsets [0] = record < 0 ? records [file_name] : record;
					
					netCDF::NcVar ncdata = files [file_name]->getVar (name.c_str ());
					if (ncdata.isNull ()) {
						ncdata = files [file_name]->addVar (name.c_str (), netcdf_type (&typeid (datatype)), scalar_dims);
						ncdata.setFill (true, NULL);
					}
					ncdata.putVar (offsets, sizes, data);
				}
		
				template <class datatype>
				static void read (const data_grid &grid, std::string file_name, std::string name, void *data, int record = -1, int flags = all_d) {
					std::vector <size_t> scalar_offsets;
					std::vector <size_t> sizes = grid.ns;
					if (flags < 0) {
						scalar_offsets = grid.offsets;
					} else {
						scalar_offsets = {grid.offsets [0]};
						sizes = {sizes [0]};
						int i = 0;
						while (i < grid.get_n_dims ()) {
							if (flags % 2 != 0) scalar_offsets.push_back (grid.offsets [i + 1]);
							if (flags % 2 != 0) sizes.push_back (grid.ns [i + 1]);
							flags = flags >> 1;
							i++;
						}
					}
					try {
						scalar_offsets [0] = record < 0 ? records [file_name] : record;
						netCDF::NcVar ncdata = files [file_name]->getVar (name.c_str ());
						if (ncdata.isNull ()) {
							throw 0;
						}
						ncdata.getVar (scalar_offsets, sizes, data);
					} catch (netCDF::exceptions::NcBadName &e) {
						failures [file_name].push_back (name);
						WARN ("Variable " << name << " not found in file");
					} catch (int &e) {
						failures [file_name].push_back (name);
						WARN ("Variable " << name << " not found in file");
					}
					
				}
			
				static bool uses_files;
			protected:
				static std::map <std::string, netCDF::NcFile *> files;
				static std::map <std::string, std::vector <netCDF::NcDim>> dims;
				static std::map <std::string, std::vector <std::string>> failures;
				static std::map <std::string, int> records;
				static std::map <std::string, bool> first;
			};
		
			/*
				TODO Accept 1D outputs
			*/
		} /* two_d */
	} /* formats */
} /* io */

#endif /* end of include guard: NETCDF_HPP_72C14A02 */
