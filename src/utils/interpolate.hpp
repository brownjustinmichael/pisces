/*!**********************************************************************
 * \file interpolate.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef INTERPOLATE_HPP_IZNK7T4T
#define INTERPOLATE_HPP_IZNK7T4T

#include "../config.hpp"

namespace utils
{
	/*!**********************************************************************
	 * \brief Interpolate over the array dy at points in
	 * 
	 * \param n The integer number of elements in the interpolation direction in dy
	 * \param m The integer number of elements perpendicular to the interpolation direction in dy
	 * \param l The integer number of elements in the interpolation direction in dx
	 * \param alpha The real multiplier on the resulting array
	 * \param x The double array of independent variables to interpolate over
	 * \param y The double array of dependent variables to interpolate over
	 * \param in The double array of independent variables to interpolate at
	 * \param out The double array of output dependent variables
	 ************************************************************************/
	template <class datatype>
	void interpolate (int n, int m, int l, datatype alpha, datatype beta, datatype *x, datatype *y, datatype *in, datatype *out, int ldy = -1, int ldout = -1) {
		TRACE ("Interpolating...");
		if (ldy == -1) {
			ldy = l;
		}
		if (ldout == -1) {
			ldout = n;
		}
		if (alpha == 0.0) {
			return;
		}
		for (int k = 0; k < n; ++k) {
			int i = 0;
			/*
				TODO Allow for reverse dx as well
			*/
			if (in [k] < x [0] || in [k] > x [l - 1]) {
				FATAL ("Interpolation out of range: " << in [k] << " not between " << x [0] << " and " << x [l - 1]);
				throw 0;
				/*
					TODO better exception?
				*/
			}
			while (in [k] > x [i]) {
				++i;
			}
			if (in [k] == x [i]) {
				for (int j = 0; j < m; ++j) {
					out [j * ldout + k] = alpha * y [j * ldy + i] + beta * out [j * ldout + k];
				}
			} else {
				for (int j = 0; j < m; ++j) {
					out [j * ldout + k] = alpha * ((y [j * ldy + i] - y [j * ldy + i - 1]) / (x [i] - x [i - 1]) * (in [k] - x [i]) + y [j * ldy + i]) + beta * out [j * ldout + k];
				}
			}
		}
	}
	
	void interpolate (int n, int m, int l, float alpha, float* x, float* y, float* in, float* out, int ldy = -1, int ldout = -1);
	
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
