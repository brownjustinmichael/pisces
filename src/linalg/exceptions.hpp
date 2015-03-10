/*!***********************************************************************
 * \file linalg/exceptions.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef EXCEPTIONS_HPP_82T7S9HZ
#define EXCEPTIONS_HPP_82T7S9HZ

#include <exception>
#include <sstream>

namespace linalg
{
	namespace exceptions
	{
		/*!**********************************************************************
		 * \brief An exception to be raised when attempting to refactorize something
		 ************************************************************************/
		class already_factorized : public std::exception
		{
		public:
			/*!**********************************************************************
			 * \brief Describe the exception
			 ************************************************************************/
			inline const char *what () const throw () {return "Matrix already factorized";}
		};
		
		/*!**********************************************************************
		 * \brief An exception to be raised when a matrix connot be factorized
		 ************************************************************************/
		class cannot_factor : public std::exception
		{
		public:
			/*!**********************************************************************
			 * \brief Describe the exception
			 ************************************************************************/
			virtual const char *what () const throw () {return "Unable to factorize matrix";}
		};
		
		/*!**********************************************************************
		 * \brief An exception to be raised when a matrix equation cannot be solved
		 ************************************************************************/
		class cannot_solve : public std::exception
		{
		public:
			/*!**********************************************************************
			 * \brief Describe the exception
			 ************************************************************************/			
			inline const char *what () const throw () {return "Unable to solve matrix";}
		};
		
		/*!**********************************************************************
		 * \brief An exception to be raised when a NaN is detected
		 ************************************************************************/
		class nan : public std::exception
		{
		public:
			/*!**********************************************************************
			 * \brief Describe the exception
			 ************************************************************************/
			inline const char *what () const throw () {return "Nan detected";}
		};
	} /* exceptions */
} /* linalg */

#endif /* end of include guard: EXCEPTIONS_HPP_82T7S9HZ */
