/*!**********************************************************************
 * \file interpolate.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef INTERPOLATE_HPP_IZNK7T4T
#define INTERPOLATE_HPP_IZNK7T4T

namespace utils
{
	/*!**********************************************************************
	 * \brief Interpolate over the array dy at x_0
	 * 
	 * \param n The integer number of elements in dx, dy
	 * \param x The double array of independent variables
	 * \param y The double array of dependent variables
	 * \param x_0 The double value at which to interpolate
	 ************************************************************************/
	void interpolate (int n, int m, double* x, double* y, double* in, double* out, int ldy = -1, int ldout = -1);
	void interpolate (int n, int m, float* x, float* y, float* in, float* out, int ldy = -1, int ldout = -1);
	
	/*!**********************************************************************
	 * \brief Interpolate a matrix multiplication at x_0
	 * 
	 * \param n The integer number of elements in x, rows in y
	 * \param x The double array of independent variables
	 * \param m The integer number of columns in y, elements in f
	 * \param y The double array matrix to interpolate
	 * \param f The double array with which to take the dot product
	 * \param x_0 The double value at which to interpolate
	 ************************************************************************/
	double dot_interpolate (int n, double* x, int m, double* y, double* f, double x_0);
	
	/*!**********************************************************************
	 * \brief Interpolate a matrix multiplication at x_0
	 * 
	 * \param n The integer number of elements in x, rows in y
	 * \param x The double array of independent variables
	 * \param m The integer number of columns in y, elements in f
	 * \param y The double array matrix to interpolate
	 * \param f The double array with which to take the dot product
	 * \param x_0 The double value at which to interpolate
	 ************************************************************************/
	float dot_interpolate (int n, float* x, int m, float* y, float* f, float x_0);
} /* utils */

#endif /* end of include guard: INTERPOLATE_HPP_IZNK7T4T */
