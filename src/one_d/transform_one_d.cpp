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
	transform <double>::transform (bases::grid <double> &i_grid, double* i_data_in, double* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags) :
	bases::plan <double> (i_element_flags, i_component_flags),
	n (i_grid.n),
	data_in (i_data_in),
	data_out (i_data_out ? i_data_out : i_data_in) {
		scalar = 1.0 / std::sqrt (2.0 * (n - 1));
		DEBUG ("FFTW Pointer: "<< data_in << " " << data_out);
		fourier_plan = fftw_plan_r2r_1d (n, data_in, data_out, FFTW_REDFT00, FFTW_ESTIMATE);
	}

	template <>
	transform <float>::transform (bases::grid <float> &i_grid, float* i_data_in, float* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags) :
	bases::plan <float> (i_element_flags, i_component_flags),
	n (i_grid.n),
	data_in (i_data_in),
	data_out (i_data_out ? i_data_out : i_data_in) {
		scalar = 1.0 / std::sqrt (2.0 * (n - 1));
		fourier_plan_float = fftwf_plan_r2r_1d (n, data_in, data_out, FFTW_REDFT00, FFTW_ESTIMATE);
	}

	template <>
	void transform <double>::execute () {
		TRACE ("Executing FFT...");

		fftw_execute (fourier_plan);
	
		for (int i = 0; i < n; ++i) {
			data_out [i] *= scalar;
		}
		
		if (!(*component_flags & transformed_vertical)) {
			*component_flags |= transformed_vertical;
		} else {
			*component_flags &= ~transformed_vertical;
		}
	}

	template <>
	void transform <float>::execute () {
		TRACE ("Executing FFT...");

		fftwf_execute (fourier_plan_float);
	
		for (int i = 0; i < n; ++i) {
			data_out [i] *= scalar;
		}
	
		if (!(*component_flags & transformed_vertical)) {
			*component_flags |= transformed_vertical;
		} else {
			*component_flags &= ~transformed_vertical;
		}
	}

	template class transform <double>;
	template class transform <float>;
} /* one_d */