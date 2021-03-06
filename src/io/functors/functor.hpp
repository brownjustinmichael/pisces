/*!**********************************************************************
 * \file functor.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FUNCTOR_HPP_DDCC8779
#define FUNCTOR_HPP_DDCC8779

#include "versions/version.hpp"

/*!**********************************************************************
 * \namespace functors
 * 
 * \brief A namespace consisting of functor objects that map one set of data to a related set
 ************************************************************************/
namespace functors
{
	/*
		TODO Make nesting functors easier
	*/
	/*!**********************************************************************
	 * \brief Abstract class for the functor object
	 * 
	 * The functor class is designed to take an instruction and apply it to some data for visualization and optimization purposes. For example, a functor could take a two dimensional grid of data and produce a one dimensional average or profile. This class serves as a wrapper for the calculate function, which returns a pointer to the processed data. Many of these functors are built with the ability to take either an array to process or another functor to allow for nested calculations (e.g. max_functor (deriv_functor) calculates the maximum derivative)
	 ************************************************************************/
	class functor
	{
	public:
		virtual ~functor () {}
		
		/*!**********************************************************************
		 * \return The version of the class
		 ************************************************************************/
		static versions::version& version () {
			static versions::version version ("1.0.1.0");
			return version;
		}
	
		/*!**********************************************************************
		 * \brief The instruction to process on the data
		 * 
		 * The class serves as a wrapper for this function
		 * 
		 * \return A pointer to the processed data, for output
		 ************************************************************************/
		virtual void *calculate () = 0;
	};
} /* functors */

#endif /* end of include guard: FUNCTOR_HPP_DDCC8779 */
