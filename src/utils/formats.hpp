/*!**********************************************************************
 * \file formats.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-06-18.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FORMATS_HPP_51078358
#define FORMATS_HPP_51078358

#include "io.hpp"

namespace io
{
	extern std::map <std::string, virtual_dump> virtual_dumps;
	
	namespace formats
	{
		namespace one_d
		{
			// /*!*******************************************************************
			//  * \brief A simple implementation of the output class
			//  * 
			//  * This class is a simple implementation of the output class.
			//  *********************************************************************/
			// class ascii : public format <ascii>
			// {
			// public:
			// 	/*!*******************************************************************
			// 	 * \param i_file_name The string name of file for output
			// 	 * \param i_n The integer number of points in the data
			// 	 * \param i_output_every An integer number of steps between outputs
			// 	 *********************************************************************/
			// 	ascii (int i_n, std::string i_comment = "#") : 
			// 	n (i_n),
			// 	comment (i_comment) {}
			// 
			// 	~ascii () {}
			// 	
			// 	std::string extension () {return ".dat";}
			// 	
			// 	virtual void open_file (std::string file_name, int file_type);
			// 
			// 	virtual void close_file ();
			// 
			// 	template <class datatype>
			// 	virtual void _write (std::string name, void *data) {
			// 		TRACE ("Writing...");
			// 		names.push_back (name);
			// 		double_data.resize (double_data.size () + 1);
			// 		float_data.resize (float_data.size () + 1);
			// 		int_data.resize (int_data.size () + 1);
			// 		if (type == &typeid (double)) {
			// 			double_data [(int) double_data.size () - 1].resize (n);
			// 			double_data [(int) double_data.size () - 1].assign ((double *) data, ((double *) data) + n);
			// 			types.push_back (&typeid (double));
			// 		} else if (type == &typeid (float)) {
			// 			float_data [(int) float_data.size () - 1].resize (n);
			// 			float_data [(int) float_data.size () - 1].assign ((float *) data, ((float *) data) + n);
			// 			types.push_back (&typeid (float));
			// 		} else if (type == &typeid (int)) {
			// 			int_data [(int) int_data.size () - 1].resize (n);
			// 			int_data [(int) int_data.size () - 1].assign ((int *) data, ((int *) data) + n);
			// 			types.push_back (&typeid (int));
			// 		}
			// 	}
			// 	virtual void write_scalar (std::string name, std::type_info *, void *) {
			// 		
			// 	}
			// 	virtual void read (std::string name, std::type_info *, void *) {
			// 		throw 0;
			// 	}
			// 	virtual void read_scalar (std::string name, std::type_info *, void *) {
			// 		throw 0;
			// 	}
			// 	
			// 	/*
			// 		TODO Write these methods
			// 	*/
			// 		
			// protected:
			// 	std::ofstream file_stream;
			// 	int n;
			// 	std::vector <const std::type_info *> types;
			// 	std::vector <std::vector <double>> double_data;
			// 	std::vector <std::vector <float>> float_data;
			// 	std::vector <std::vector <int>> int_data;
			// 	std::vector <std::string> names;
			// 	std::string comment;
			// };
		
			/*
				TODO Make legible headers
			*/
		
			// class netcdf : public format
			// {
			// public:
			// 	netcdf (int i_n, int i_n_offset = 0) :
			// 	n (i_n), n_offset (i_n_offset) {}
			// 
			// 	virtual ~netcdf () {}
			// 	
			// 	std::string extension () {return ".cdf";}
			// 	
			// 	virtual void open_file (std::string file_name, int file_type);
			// 
			// 	virtual void close_file ();
			// 
			// 	virtual void write (std::string name, double *);
			// 	virtual void write (std::string name, float *);
			// 	virtual void write (std::string name, int *);
			// 	
			// 	virtual void write_scalar (std::string name, double *);
			// 	virtual void write_scalar (std::string name, float *);
			// 	virtual void write_scalar (std::string name, int *);
			// 	
			// 	virtual void read (std::string name, double *);
			// 	virtual void read (std::string name, float *);
			// 	virtual void read (std::string name, int *);
			// 	
			// 	virtual void read_scalar (std::string name, double *);
			// 	virtual void read_scalar (std::string name, float *);
			// 	virtual void read_scalar (std::string name, int *);
			// 
			// protected:
			// 	int n;
			// 	int n_offset;
			// 	void * datafile;
			// 	netCDF::NcDim zDim;
			// };
		} /* one_d */
	
		namespace two_d
		{
			class virtual_format
			{
			public:
				virtual_format () {}
		
				~virtual_format () {}
		
				static std::string extension () {return "";}
			
				static void open_file (std::string file_name, int file_type, int n_max, int m_max, int l_max) {
					DEBUG ("Opening virtual file.");
					// if (file_type == read_file && (!virtual_dumps [file_name])) {
						// ERROR ("Virtual file doesn't exist.");
						// throw 0;
					// } 
					/*
						TODO Check for virtual file existence
					*/
					virtual_dumps [file_name];
				}
			
				static void close_file (std::string file_name, int file_type) {
					DEBUG ("Closing virtual file.");
				}
			
				template <class datatype>
				static void write (std::string file_name, std::string name, datatype *data, int n = 1, int m = 1, int l = 1, int n_offset = 0, int m_offset = 0, int l_offset = 0, int record = -1) {
					DEBUG ("Writing to virtual file");
					std::stringstream debug;
					virtual_dumps [file_name].add_var <datatype> (name, n, m);
					virtual_dumps [file_name].put <datatype> (name, (datatype *) data, n, m);
					DEBUG ("Dimensions " << n << ", " << m << ", " << l);
					for (int i = 0; i < n; ++i) {
						for (int j = 0; j < m; ++j) {
							debug << data [i * m + j] << ":" << virtual_dumps [file_name].index <datatype> (name, i, j) << " ";
						}
						if (file_name == "main/dump") DEBUG (debug.str ());
						debug.str ("");
					}
				}
		
				template <class datatype>
				static void write_scalar (std::string file_name, std::string name, datatype *data, int record = -1) {
					virtual_dumps [file_name].add_var <datatype> (name);
					virtual_dumps [file_name].put <datatype> (name, (datatype *) data);
				}
		
				template <class datatype>
				static void read (std::string file_name, std::string name, datatype *data, int n = 1, int m = 1, int l = 1, int n_offset = 0, int m_offset = 0, int l_offset = 0, int record = -1) {
					DEBUG ("Reading from virtual file");
					virtual_dumps [file_name].get <datatype> (name, (datatype *) data, n, m);
					for (int i = 0; i < n; ++i) {
						for (int j = 0; j < m; ++j) {
							if (*(data + i * m + j) != *(data + i * m + j)) {
								FATAL ("NaN read in." << i << " " << j << " " << *(data + i * m + j));
								throw 0;
							}
						}
					}
				}
		
				template <class datatype>
				static void read_scalar (std::string file_name, std::string name, datatype *data, int record = -1) {
					virtual_dumps [file_name].get <datatype> (name, (datatype *) data);
				}
			};
			
		
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
						DEBUG ("Read " << *data);
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

#endif /* end of include guard: FORMATS_HPP_51078358 */
