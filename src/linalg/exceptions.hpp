/*!***********************************************************************
 * \file exceptions.hpp
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
	} /* exceptions */
} /* linalg */

#endif /* end of include guard: EXCEPTIONS_HPP_82T7S9HZ */
