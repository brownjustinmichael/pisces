/*!**********************************************************************
 * \file format.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FORMAT_HPP_6ADF5F85
#define FORMAT_HPP_6ADF5F85

namespace io
{
	/*!**********************************************************************
	 * \brief A set of io flags to be used with the input/output classes specifying the file type
	 ************************************************************************/
	enum io_flags {
		read_file = 0,
		replace_file = 1,
		append_file = 2
	};
} /* io */

#endif /* end of include guard: FORMAT_HPP_6ADF5F85 */
