/*!**********************************************************************
 * \file interpolate.hpp
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef INTERPOLATE_HPP_IZNK7T4T
#define INTERPOLATE_HPP_IZNK7T4T

namespace utils
{
	double interpolate (int n, double* dx, double* dy, double x);
	
	void matrix_interpolate (int n, double* dx, int m, double* dy, double da, int incy, double* douty, double x);
	
} /* utils */

#endif /* end of include guard: INTERPOLATE_HPP_IZNK7T4T */
