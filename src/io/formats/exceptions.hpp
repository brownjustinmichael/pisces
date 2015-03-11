/*!***********************************************************************
 * \file io/formats/exceptions.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef EXCEPTIONS_HPP_82T7S9HQ
#define EXCEPTIONS_HPP_82T7S9HQ

#include <exception>
#include <sstream>

/*!**********************************************************************
 * \namespace formats::exceptions
 * 
 * \brief A namespace containing the exceptions for the format classes
 ************************************************************************/
namespace formats
{
	namespace exceptions
	{
		/*!*******************************************************************
		 * \brief This is an exception that can occur when opening files
		 *********************************************************************/
		class file_exception : public std::exception
		{
		private:
			std::string file_name; //!< The file name that has the issue
			
		public:
			/*!**********************************************************************
			 * \param i_file_name The file name that has the issue
			 ************************************************************************/
			file_exception (std::string i_file_name) {
				file_name = i_file_name;
			}
			
			~file_exception () throw () {}
			
			/*!*******************************************************************
			 * \brief Describe the nature of the exception
			 *********************************************************************/
			inline const char *what () const throw () {
				std::stringstream message;
				message << "File " << file_name << " couldn't be opened";
				return message.str ().c_str ();
			}
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
		
		/*!**********************************************************************
		 * \brief This exception occurs when a variable is missing from a file
		 ************************************************************************/
		class bad_variables : public std::exception
		{
		public:
			/*!**********************************************************************
			 * \brief Describe the nature of the exception
			 ************************************************************************/
			inline const char *what () const throw () {return "Bad variable(s)";}
		};
	} /* exceptions */
} /* formats */

#endif /* end of include guard: EXCEPTIONS_HPP_82T7S9HQ */
