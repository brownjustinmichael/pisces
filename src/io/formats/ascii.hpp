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
#include <memory>
#include <string>
#include <fstream>

#include "logger/logger.hpp"

#include "format.hpp"
#include "exceptions.hpp"

namespace io
{
	namespace formats
	{
		/*!**********************************************************************
		 * \brief A file format designed to interface with ASCII files
		 * 
		 * \copydoc virtual_format
		 ************************************************************************/
		class ascii
		{
		private:
			/*
				TODO Different files should have different comment strings
			*/
			static std::string comment; //!< The comment string for the file
			static std::map <std::string, std::shared_ptr <std::ofstream>> file_streams; //!< A map of streams to convert from file name to file stream
			static std::map <std::string, int> count; //!< A map of ints to count the number of times a file has been opened
			static std::map <std::string, int> file_types; //!< A map of the io flags for each file
			static std::map <std::string, std::shared_ptr <std::stringstream>> header; //!< The header for each file
			static std::map <std::string, std::shared_ptr <std::stringstream>> body; //!< The current line of each file
			
		public:
			static bool print_headers; //!< A boolean determining whether or not to print the headers of a file
			static bool uses_files; //!< A boolean determining whether a format uses files (true)
			
			ascii () {}
			
			~ascii () {}
			
			/*!**********************************************************************
			 * \copydoc virtual_format::extension
			 ************************************************************************/
			static std::string extension () {return ".dat";}
			
			/*!**********************************************************************
			 * \copydoc virtual_format::open_file
			 ************************************************************************/
			static void open_file (const data_grid &grid, std::string file_name, int file_type) {
				TRACE ("Opening ASCII file");
				if (!(file_streams [file_name])) {
					file_streams [file_name] = std::shared_ptr <std::ofstream> (new std::ofstream);
				}
				if (file_type == read_file) {
					file_streams [file_name]->open (file_name, std::ostream::in);
				} else if (file_type == replace_file || count [file_name] == 0) {
					file_streams [file_name]->open (file_name, std::ostream::out | std::ostream::trunc);
				} else if (file_type == append_file) {
					if (!(file_streams [file_name]->is_open ())) {
						file_streams [file_name]->open (file_name, std::ostream::out | std::ostream::app);
					}
				} else {
					ERROR ("Unknown file type");
					throw 0;
				}
				
				DEBUG ("HERE");
				
				if (!file_streams [file_name]->is_open ()) {
					ERROR ("Failed to open file " << file_name);
					throw exceptions::file_exception (file_name);
				}
				
				if (!header [file_name]) {
					header [file_name] = std::shared_ptr <std::stringstream> (new std::stringstream);
				}
				// if (header [file_name]->str () != "") {
				// 	*file_streams [file_name] << comment << " " << header [file_name]->str () << "\n";
				// }
				if (!body [file_name]) {
					body [file_name] = std::shared_ptr <std::stringstream> (new std::stringstream);
				}
				// if (body [file_name]->str () != "") {
				// 	*file_streams [file_name] << body [file_name]->str () << "\n";
				// }
				
				count [file_name] += 1;
				file_types [file_name] = file_type;
				
				header [file_name]->str ("");
				body [file_name]->str ("");
				
				DEBUG ("DONE");
			}
			
			/*!**********************************************************************
			 * \copydoc virtual_format::close_file
			 ************************************************************************/
			static void close_file (std::string file_name, int file_type) {
				TRACE ("Closing file");
				if (print_headers && header [file_name]->str () != "") {
					*file_streams [file_name] << comment << " " << header [file_name]->str () << "\n";
				}
				*file_streams [file_name] << body [file_name]->str () << "\n";
				file_streams [file_name]->close ();
				TRACE ("Done");
			}
			
			/*!**********************************************************************
			 * \copydoc virtual_format::write
			 ************************************************************************/
			template <class datatype>
			static void write (const data_grid &grid, std::string file_name, std::string name, void *data, int record = -1, int flags = all_d) {
				if (flags != scalar) {
					FATAL ("ASCII write not implemented");
					throw 0;
				}
				if (file_types [file_name] != append_file || count [file_name] == 1) {
					*header [file_name] << name << "\t";
				}
				*body [file_name] << * (datatype *) data << "\t";
			}
			
			/*!**********************************************************************
			 * \copydoc virtual_format::read
			 ************************************************************************/
			template <class datatype>
			static void read (const data_grid &grid, std::string file_name, std::string name, void *data, int record = -1, int flags = all_d) {
				FATAL ("ASCII read not implemented");
				throw 0;
			}
		};
	} /* formats */
} /* io */

#endif /* end of include guard: FORMATS_HPP_51078358 */
