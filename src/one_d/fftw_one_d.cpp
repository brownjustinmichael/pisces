/*!**********************************************************************
 * \file fftw_one_d.cpp
 * /Users/justinbrown/Dropbox/spectral_element/src
 * 
 * Created by Justin Brown on 2013-08-20.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "fftw_one_d.hpp"
#include <cmath>

namespace one_d
{
	template <>
	fftw_cosine <double>::fftw_cosine (bases::element <double>* i_element_ptr, int i_n, int i_name_in, int i_name_out, int i_flags) : 
	bases::explicit_plan <double> (i_element_ptr, i_n, i_name_in, i_name_out, i_flags) {
		scalar = 1.0 / std::sqrt (2.0 * (n - 1));
		fourier_plan = fftw_plan_r2r_1d (n, data_in, data_out, FFTW_REDFT00, FFTW_ESTIMATE);
	}
	
	template <>
	fftw_cosine <float>::fftw_cosine (bases::element <float>* i_element_ptr, int i_n, int i_name_in, int i_name_out, int i_flags) : 
	bases::explicit_plan <float> (i_element_ptr, i_n, i_name_in, i_name_out, i_flags) {
		scalar = 1.0 / std::sqrt (2.0 * (n - 1));
		fourier_plan_float = fftwf_plan_r2r_1d (n, data_in, data_out, FFTW_REDFT00, FFTW_ESTIMATE);
	}
	
	template <>
	void fftw_cosine <double>::execute () {
		TRACE ("Executing...");
		
		bases::explicit_plan <double>::execute ();
		
		fftw_execute (fourier_plan);
		
		for (int i = 0; i < n; ++i) {
			data_out [i] *= scalar;
		}
		
	}
	
	template <>
	void fftw_cosine <float>::execute () {
		TRACE ("Executing...");
		
		bases::explicit_plan <float>::execute ();
		
		fftwf_execute (fourier_plan_float);
		
		for (int i = 0; i < n; ++i) {
			data_out [i] *= scalar;
		}
		
	}
	
	template class fftw_cosine <double>;
	template class fftw_cosine <float>;
} /* one_d */