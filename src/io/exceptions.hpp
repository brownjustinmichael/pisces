// 
//! \file exceptions.hpp
//  io
//  
//  Created by Justin Brown on 2013-04-01.
//  Copyright 2013 Justin Brown. All rights reserved.
// 

#ifndef EXCEPTIONS_HPP_82T7S9HQ
#define EXCEPTIONS_HPP_82T7S9HQ

#include <exception>

namespace io
{
	namespace exceptions
	{
		class file_exception : public std::exception
		{
		public:
			inline const char *what () {return "File couldn't be opened";}
		};
	} /* exceptions */
} /* io */

#endif /* end of include guard: EXCEPTIONS_HPP_82T7S9HQ */
