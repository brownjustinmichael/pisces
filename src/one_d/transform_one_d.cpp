/*!**********************************************************************
 * \file transform_one_d.cpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-20.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "transform_one_d.hpp"
#include <cmath>

namespace one_d
{
	template <>
	fftw_cosine <double>::fftw_cosine (bases::grid <double> &i_grid, double* i_data_in, double* i_data_out) :
	n (i_grid.n),
	data_in (i_data_in),
	data_out (i_data_out ? i_data_out : i_data_in) {
		scalar = 1.0 / std::sqrt (2.0 * (n - 1));
		fourier_plan = fftw_plan_r2r_1d (n, data_in, data_out, FFTW_REDFT00, FFTW_ESTIMATE);
	}
	
	template <>
	fftw_cosine <float>::fftw_cosine (bases::grid <float> &i_grid, float* i_data_in, float* i_data_out) :
	n (i_grid.n),
	data_in (i_data_in),
	data_out (i_data_out ? i_data_out : i_data_in) {
		scalar = 1.0 / std::sqrt (2.0 * (n - 1));
		fourier_plan_float = fftwf_plan_r2r_1d (n, data_in, data_out, FFTW_REDFT00, FFTW_ESTIMATE);
	}
	
	template <>
	void fftw_cosine <double>::execute (int &element_flags, int &component_flags) {
		TRACE ("Executing FFT...");

		fftw_execute (fourier_plan);
		
		for (int i = 0; i < n; ++i) {
			data_out [i] *= scalar;
		}
		
	}
	
	template <>
	void fftw_cosine <float>::execute (int &element_flags, int &component_flags) {
		TRACE ("Executing FFT...");
		
		fftwf_execute (fourier_plan_float);
		
		for (int i = 0; i < n; ++i) {
			data_out [i] *= scalar;
		}
		
	}
	
	template class fftw_cosine <double>;
	template class fftw_cosine <float>;
} /* one_d */