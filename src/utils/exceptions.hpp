/*!***********************************************************************
 * \file exceptions.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef EXCEPTIONS_HPP_82T7S9HQ
#define EXCEPTIONS_HPP_82T7S9HQ

#include <exception>

namespace io
{
	namespace exceptions
	{
		/*!*******************************************************************
		 * \brief This is an exception that can occur when opening files
		 *********************************************************************/
		class file_exception : public std::exception
		{
		public:
			/*!*******************************************************************
			 * \brief Describe the nature of the exception
			 *********************************************************************/
			inline const char *what () {return "File couldn't be opened";}
		};
	} /* exceptions */
} /* io */

#endif /* end of include guard: EXCEPTIONS_HPP_82T7S9HQ */
