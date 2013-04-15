/*!***********************************************************************
 * \file fft.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-15.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FFT_HPP_P3TP70YE
#define FFT_HPP_P3TP70YE

#include <fftw3.h>
#include "../config.hpp"
#include "../plan.hpp"

namespace fft
{
	class fft : public plan
	{
	public:
		virtual ~fft () {}
		virtual void execute (double timestep, int *flags) = 0;
	};
	
	class fftw_cosine : public fft
	{
	public:
		fftw_cosine (int i_n, double *i_data_in, double *i_data_out, int i_fftw_flags = FFTW_ESTIMATE) {
			n = i_n,
			data_out = i_data_out;
			scalar = 1.0 / sqrt (2.0 * (i_n - 1));
			
			fourier_plan = fftw_plan_r2r_1d (i_n + 1, i_data_in, i_data_out, FFTW_REDFT00, i_fftw_flags);
		}
		virtual ~fftw_cosine () {}
		
		void execute (double timestep, int *flags) {
			int i;
			
			fftw_execute (fourier_plan);
			
			for (i = 0; i < n; ++i) {
				data_out [i] *= scalar;				
			}
		}
	
	private:
		int n;
		double *data_out;
		
		double scalar;
		fftw_plan fourier_plan;
	};
} /* fft */

#endif /* end of include guard: FFT_HPP_P3TP70YE */
