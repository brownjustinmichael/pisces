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
	fftw_cosine <double>::fftw_cosine (int i_n, double* i_data_in, double* i_data_out) : 
	explicit_plan <double> (i_n, i_data_in, i_data_out) {
		scalar = 1.0 / std::sqrt (2.0 * (n - 1));
		fourier_plan = fftw_plan_r2r_1d (n, data_in, data_out, FFTW_REDFT00, FFTW_ESTIMATE);
	}
	
	template <>
	fftw_cosine <float>::fftw_cosine (int i_n, float* i_data_in, float* i_data_out) : 
	explicit_plan <float> (i_n, i_data_in, i_data_out) {
		scalar = 1.0 / std::sqrt (2.0 * (n - 1));
		fourier_plan_float = fftwf_plan_r2r_1d (n, data_in, data_out, FFTW_REDFT00, FFTW_ESTIMATE);
	}
	
	template <>
	void fftw_cosine <double>::execute () {
		TRACE ("Executing...");
		
		fftw_execute (fourier_plan);
		
		for (int i = 0; i < n; ++i) {
			data_out [i] *= scalar;
		}
		
	}
	
	template <>
	void fftw_cosine <float>::execute () {
		TRACE ("Executing...");
		
		fftwf_execute (fourier_plan_float);
		
		for (int i = 0; i < n; ++i) {
			data_out [i] *= scalar;
		}
		
	}
	
	template class fftw_cosine <double>;
	template class fftw_cosine <float>;
} /* one_d */