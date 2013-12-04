/*!**********************************************************************
 * \file interpolate.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "interpolate.hpp"
#include "../config.hpp"

namespace utils
{
	void interpolate (int n, int m, int l, double alpha, double* x, double* y, double* in, double* out, int ldy, int ldout) {
		TRACE ("Interpolating...");
		if (ldy == -1) {
			ldy = m;
		}
		if (ldout == -1) {
			ldout = n;
		}
		for (int k = 0; k < n; ++k) {
			int i = 1;
			/*
				TODO Allow for reverse dx as well
			*/
			if (in [k] < x [0] || in [k] > x [l - 1]) {
				FATAL ("Interpolation out of range.");
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
					out [j * ldout + k] += alpha * y [j * ldy + i];
				}
			} else {
				for (int j = 0; j < m; ++j) {
					out [j * ldout + k] += alpha * ((y [j * ldy + i] - y [j * ldy + i - 1]) / (x [i] - x [i - 1]) * (in [k] - x [i]) + y [j * ldy + i]);
				}
			}
		}
	}
	
	void interpolate (int n, int m, int l, float alpha, float* x, float* y, float* in, float* out, int ldy, int ldout) {
		TRACE ("Interpolating...");
		if (ldy == -1) {
			ldy = m;
		}
		if (ldout == -1) {
			ldout = n;
		}
		for (int k = 0; k < n; ++k) {
			int i = 1;
			/*
				TODO Allow for reverse dx as well
			*/
			if (in [k] < x [0] || in [k] > x [l - 1]) {
				FATAL ("Interpolation out of range.");
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
					out [j * ldy + k] += alpha * y [j * ldy + i];
				}
			} else {
				for (int j = 0; j < m; ++j) {
					out [j * ldy + k] += alpha * ((y [j * ldy + i] - y [j * ldy + i - 1]) / (x [i] - x [i - 1]) * (in [j] - x [i]) + y [j * ldy + i]);
				}
			}
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