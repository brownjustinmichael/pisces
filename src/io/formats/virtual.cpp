/*!**********************************************************************
 * \file virtual.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "virtual.hpp"

namespace formats
{
	std::map <std::string, formats::virtual_file> virtual_files;
	
	bool virtual_format::uses_files = false;
} /* formats */
