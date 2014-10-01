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
			
				static void open_file (std::string file_name, int file_type, int n_max, int m_max, int l_max);
		
				static void close_file (std::string file_name, int file_type);
			
				template <class datatype>
				static void write (std::string file_name, std::string name, datatype *data, int n = 1, int m = 1, int l = 1, int n_offset = 0, int m_offset = 0, int l_offset = 0, int record = -1) {
					std::vector <size_t> offsets = {(size_t) (record < 0 ? records [file_name] : record), (size_t) n_offset, (size_t) m_offset};
					std::vector <size_t> sizes = {(size_t) 1, (size_t) n, (size_t) m};
					netCDF::NcVar ncdata = files [file_name]->getVar (name.c_str ());
					if (ncdata.isNull ()) {
						ncdata = files [file_name]->addVar (name.c_str (), netcdf_type (&typeid (datatype)), dims [file_name]);
						ncdata.setFill (true, NULL);
					}
					ncdata.putVar (offsets, sizes, data);
				}
		
				template <class datatype>
				static void write_scalar (std::string file_name, std::string name, datatype *data, int n = 1, int m = 1, int l = 1, int n_offset = 0, int m_offset = 0, int l_offset = 0, int record = -1) {
					std::vector <netCDF::NcDim> scalar_dims = {dims [file_name] [0]};
					// std::vector <netCDF::NcDim> scalar_dims;
					std::vector <size_t> scalar_offset = {(size_t) (record < 0 ? records [file_name] : record)};
					// std::vector <size_t> scalar_offset;
					std::vector <size_t> sizes = {(size_t) 1};
					netCDF::NcVar ncdata = files [file_name]->getVar (name.c_str ());
					if (ncdata.isNull ()) {
						ncdata = files [file_name]->addVar (name.c_str (), netcdf_type (&typeid (datatype)), scalar_dims);
					}
					ncdata.putVar (scalar_offset, sizes, data);
				}
		
				template <class datatype>
				static void read (std::string file_name, std::string name, datatype *data, int n = 1, int m = 1, int l = 1, int n_offset = 0, int m_offset = 0, int l_offset = 0, int record = -1) {
					try {
						std::vector <size_t> offsets = {(size_t) (record < 0 ? records [file_name] : record), (size_t) n_offset, (size_t) m_offset};
						std::vector <size_t> sizes = {(size_t) 1, (size_t) n, (size_t) m};
						netCDF::NcVar ncdata = files [file_name]->getVar (name.c_str ());
						if (ncdata.isNull ()) {
							throw 0;
						}
						ncdata.getVar (offsets, sizes, data);
						for (int i = 0; i < n; ++i) {
							for (int j = 0; j < m; ++j) {
								if (*(data + i * m + j) != *(data + i * m + j)) {
									FATAL ("NaN read in. ");
									throw 0;
								}
							}
						}
					} catch (netCDF::exceptions::NcBadName &e) {
						failures [file_name].push_back (name);
						WARN ("Variable " << name << " not found in file");
					} catch (int &e) {
						failures [file_name].push_back (name);
						WARN ("Variable " << name << " not found in file");
					}
				}
		
				template <class datatype>
				static void read_scalar (std::string file_name, std::string name, datatype *data, int record = -1) {
					try {
						std::vector <size_t> scalar_offset = {(size_t) (record < 0 ? records [file_name] : record), (size_t) 0,  (size_t) 0};
						std::vector <size_t> sizes = {(size_t) 1, (size_t) 0, (size_t) 0};
						netCDF::NcVar ncdata = files [file_name]->getVar (name.c_str ());
						if (ncdata.isNull ()) {
							throw 0;
						}
						ncdata.getVar (scalar_offset, sizes, data);
					} catch (netCDF::exceptions::NcBadName &e) {
						failures [file_name].push_back (name);
						WARN ("Variable " << name << " not found in file");
					} catch (int &e) {
						failures [file_name].push_back (name);
						WARN ("Variable " << name << " not found in file");
					}
				}
			
			protected:
				static std::map <std::string, netCDF::NcFile *> files;
				static std::map <std::string, std::vector <netCDF::NcDim>> dims;
				static std::map <std::string, std::vector <std::string>> failures;
				static std::map <std::string, int> records;
			};
		
			/*
				TODO Accept 1D outputs
			*/
		} /* two_d */
	} /* formats */
} /* io */

#endif /* end of include guard: NETCDF_HPP_72C14A02 */
