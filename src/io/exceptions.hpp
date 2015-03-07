/*!***********************************************************************
 * \file exceptions.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef EXCEPTIONS_HPP_82T7S9HJ
#define EXCEPTIONS_HPP_82T7S9HJ

#include <exception>
#include <sstream>
#include <iostream>

namespace io
{
	namespace exceptions
	{
		/*!**********************************************************************
		 * \brief An exception to be raised when a parameter key does not exist
		 ************************************************************************/
		class key_does_not_exist : public std::exception
		{
		public:
			std::string key_name; //!< The key that does not exist
			
			/*!**********************************************************************
			 * \param i_key_name The key that does not exist
			 ************************************************************************/
			key_does_not_exist (std::string i_key_name) :
			key_name (i_key_name) {
				std::cout << i_key_name;
			}
			
			~key_does_not_exist () throw () {}
			
			/*!**********************************************************************
			 * \brief Explain the exception
			 ************************************************************************/
			inline const char *what () const throw () {
				std::stringstream message;
				message << "Key " << key_name << " not found";
				return message.str ().c_str ();
			}
		};
	} /* exceptions */
} /* io */

#endif /* end of include guard: EXCEPTIONS_HPP_82T7S9HJ */
