/*!**********************************************************************
 * \file functor.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FUNCTOR_HPP_DDCC8779
#define FUNCTOR_HPP_DDCC8779

namespace io
{
	namespace functors
	{
		/*!**********************************************************************
		 * \brief Abstract class for the functor object
		 * 
		 * The functor class is designed to take an instruction and apply it to some data for visualization and optimization purposes. For example, a functor could take a two dimensional grid of data and produce a one dimensional average or profile. This class serves as a wrapper for the calculate function, which returns a pointer to the processed data.
		 ************************************************************************/
		template <class datatype>
		class functor
		{
		public:
			virtual ~functor () {}
		
			/*!**********************************************************************
			 * \brief The instruction to process on the data
			 * 
			 * The class serves as a wrapper for this function
			 * 
			 * \return A pointer to the processed data, for output
			 ************************************************************************/
			virtual datatype *calculate () = 0;
		};
	} /* functors */
} /* io */

#endif /* end of include guard: FUNCTOR_HPP_DDCC8779 */
