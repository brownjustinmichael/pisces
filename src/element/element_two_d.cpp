/*!**********************************************************************
 * \file element_two_d.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-12.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "element_two_d.hpp"

namespace two_d
{
	namespace fourier
	{
		namespace cosine
		{
			template <class datatype>
			int element <datatype>::mode = mode_flag;
			
			template class element <double>;
		} /* cosine */
		
		namespace chebyshev
		{
			template <class datatype>
			int element <datatype>::mode = mode_flag;

			template class element <double>;
		} /* chebyshev */
	} /* fourier */
} /* two_d */