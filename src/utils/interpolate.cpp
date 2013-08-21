/*!**********************************************************************
 * \file interpolate.cpp
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "interpolate.hpp"
#include "../config.hpp"

namespace utils
{
	double interpolate (int n, double* dx, double* dy, double x) {
		TRACE ("Interpolating...");
		int i = 1;
		/*
			TODO Allow for reverse dx as well
		*/
		if (x < dx [0] || x > dx [n - 1]) {
			throw 0;
			/*
				TODO better exception?
			*/
		}
		while (x < dx [i]) {
			++i;
		}
		if (x == dx [i]) {
			return dy [i];
		} else {
			return (dy [i] - dy [i - 1]) / (dx [i] - dx [i - 1]) * (x - dx [i]) + dy [i];
		}
	}
	
	double dot_interpolate (int n, double* dx, int m, double* dy, double* df, double x) {
		TRACE ("Interpolating matrices...");
		int i = 1;
		double running_sum = 0.0;
		/*
			TODO Allow for reverse dx as well
		*/
		if (x < dx [0] || x > dx [n - 1]) {
			throw 0;
			/*
				TODO better exception?
			*/
		}
		while (x > dx [i]) {
			++i;
		}
		if (x == dx [i]) {
			for (int j = 0; j < m; ++j) {
				running_sum += df [j] * dy [j * m + i];
			}
		} else {
			for (int j = 0; j < m; ++j) {
				running_sum += df [j] * ((dy [j * m + i] - dy [j * m + i - 1]) / (dx [i] - dx [i - 1]) * (x - dx [i]) + dy [j * m + i]);
			}
		}
		return running_sum;
	}
	
	float dot_interpolate (int n, float* dx, int m, float* dy, float* df, float x) {
		TRACE ("Interpolating matrices...");
		int i = 1;
		float running_sum = 0.0;
		/*
			TODO Allow for reverse dx as well
		*/
		if (x < dx [0] || x > dx [n - 1]) {
			throw 0;
			/*
				TODO better exception?
			*/
		}
		while (x > dx [i]) {
			++i;
		}
		if (x == dx [i]) {
			for (int j = 0; j < m; ++j) {
				running_sum += df [j] * dy [j * m + i];
			}
		} else {
			for (int j = 0; j < m; ++j) {
				running_sum += df [j] * ((dy [j * m + i] - dy [j * m + i - 1]) / (dx [i] - dx [i - 1]) * (x - dx [i]) + dy [j * m + i]);
			}
		}
		return running_sum;
	}
	
} /* utils */