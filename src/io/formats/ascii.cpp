/*!**********************************************************************
 * \file formats.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-06-18.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "ascii.hpp"

namespace io
{
	namespace formats
	{
		std::string ascii::comment = "#";
		bool ascii::print_headers = true;
		std::map <std::string, std::shared_ptr <std::ofstream>> ascii::file_streams;
		std::map <std::string, int> ascii::count;
		std::map <std::string, int> ascii::file_types;
		std::map <std::string, std::shared_ptr <std::stringstream>> ascii::header;
		std::map <std::string, std::shared_ptr <std::stringstream>> ascii::body;
		bool ascii::uses_files = true;
	} /* formats */
} /* io */