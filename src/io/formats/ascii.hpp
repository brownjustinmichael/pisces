/*!**********************************************************************
 * \file formats.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-06-18.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FORMATS_HPP_51078358
#define FORMATS_HPP_51078358

#include <map>
#include <string>
#include <fstream>

#include "logger/logger.hpp"

#include "format.hpp"
#include "exceptions.hpp"

namespace io
{
	namespace formats
	{
		class ascii
		{
		private:
			static std::string comment;
			static std::map <std::string, std::ofstream> file_streams;
			static std::map <std::string, int> count;
			static std::map <std::string, int> file_types;
			static std::map <std::string, std::stringstream> header;
			static std::map <std::string, std::stringstream> body;
			
		public:
			ascii () {}
			
			~ascii () {}
			
			static std::string extension () {return ".dat";}
			
			static void open_file (const data_grid &grid, std::string file_name, int file_type) {
				if (file_type == read_file) {
					file_streams [file_name].open (file_name, std::ostream::in);
				} else if (file_type == replace_file) {
					file_streams [file_name].open (file_name, std::ostream::out | std::ostream::trunc);
				} else if (file_type == append_file) {
					if (!(file_streams [file_name].is_open ())) {
						file_streams [file_name].open (file_name, std::ostream::out | std::ostream::trunc);
					}
				} else {
					ERROR ("Unknown file type");
					throw 0;
				}
				
				if (! file_streams [file_name].is_open ()) {
					ERROR ("Failed to open file " << file_name);
					throw exceptions::file_exception ();
				}
				
				if (header [file_name].str () != "") {
					file_streams [file_name] << comment << " " << header [file_name].str () << "\n";
				}
				if (body [file_name].str () != "") {
					file_streams [file_name] << body [file_name].str () << "\n";
				}
				
				count [file_name] += 1;
				file_types [file_name] = file_type;
				
				header [file_name].str ("");
				body [file_name].str ("");
			}
			
			static void close_file (std::string file_name, int file_type) {
				if (header [file_name].str () != "") {
					file_streams [file_name] << comment << " " << header [file_name].str () << "\n";
				}
				file_streams [file_name] << body [file_name].str () << "\n";
				file_streams [file_name].close ();
			}
			
			template <class datatype>
			static void write (const data_grid &grid, std::string file_name, std::string name, void *data, int record = -1, int flags = all_d) {
				if (flags != scalar) {
					FATAL ("ASCII write not implemented");
					throw 0;
				}
				if (file_types [file_name] != append_file || count [file_name] == 1) {
					header [file_name] << name << " ";
				}
				body [file_name] << * (datatype *) data << " ";
			}
			
			template <class datatype>
			static void read (const data_grid &grid, std::string file_name, std::string name, void *data, int record = -1, int flags = all_d) {
				FATAL ("ASCII read not implemented");
				throw 0;
			}
			
			template <class datatype>
			static void read_scalar (std::string file_name, std::string name, datatype *data, int record = -1) {
				FATAL ("ASCII read scalar not implemented");
				throw 0;
			}
		};
	} /* formats */
} /* io */

#endif /* end of include guard: FORMATS_HPP_51078358 */