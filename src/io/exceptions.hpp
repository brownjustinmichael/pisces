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
		inline const char *what () const throw () {return "File couldn't be opened";}
	};

	/*!**********************************************************************
	 * \brief This exception occurs when a function is called that can't handle the given type
	 ************************************************************************/
	class bad_type : public std::exception
	{
	public:
		/*!**********************************************************************
		 * \brief Describe the nature of the exception
		 ************************************************************************/
		inline const char *what () const throw () {return "Bad type passed";}
	};

	class key_does_not_exist : public std::exception
	{
	public:
		key_does_not_exist (std::string i_key_name) :
		key_name (i_key_name) {}

		~key_does_not_exist () throw () {}

		inline const char *what () const throw () {
			std::stringstream message;
			message << "Key " << key_name << " not found";
			return message.str ().c_str ();
		}

	// private:
		std::string key_name;
	};
	
	namespace io
	{
		class bad_variables : public std::exception
		{
		public:
			inline const char *what () const throw () {return "Bad variable(s)";}
		};
	} /* io */
} /* exceptions */

#endif /* end of include guard: EXCEPTIONS_HPP_82T7S9HJ */
