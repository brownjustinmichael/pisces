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
	
	class already_factorized : public std::exception
	{
	public:
		inline const char *what () const throw () {return "Matrix already factorized";}
	};
	
	class cannot_factor : public std::exception
	{
	public:
		virtual const char *what () const throw () {return "Unable to factor matrix";}
	};
	
	class cannot_solve : public std::exception
	{
	public:
		inline const char *what () const throw () {return "Unable to solve matrix";}
	};
	
	class nan : public std::exception
	{
	public:
		inline const char *what () const throw () {return "Nan detected";}
	};
	
	class mpi_fatal : public std::exception
	{
	public:
		inline const char *what () const throw () {return "One or more MPI Processes failed";}
	};
	
	class mesh_adapt : public std::exception
	{
	public:
		inline const char *what () const throw () {return "Mesh adaptation needed";}
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

#endif /* end of include guard: EXCEPTIONS_HPP_82T7S9HQ */
