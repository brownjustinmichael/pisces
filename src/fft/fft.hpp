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
	class fft : public explicit_plan
	{
	public:
		fft (int i_n, double *i_data_in, double *i_data_out, int *i_flags_ptr = NULL) : explicit_plan (i_n, i_data_in, i_data_out, i_flags_ptr) {}
		virtual ~fft () {}
		virtual void execute () = 0;
	};
	
	class fftw_cosine : public fft
	{
	public:
		fftw_cosine (int i_n, double *i_data_in, double *i_data_out, int *i_flags_ptr = NULL) : fft (i_n, i_data_in, i_data_out, i_flags_ptr) {
			TRACE ("Instantiating...");
			
			INFO ("FFTW_ESTIMATE = " << FFTW_ESTIMATE);
			
			scalar = 1.0 / sqrt (2.0 * (i_n - 1));
			
			fourier_plan = fftw_plan_r2r_1d (i_n + 1, i_data_in, i_data_out, FFTW_REDFT00, *flags_ptr);
			
			TRACE ("Instantiated.")
		}
		virtual ~fftw_cosine () {}
		
		void execute () {
			TRACE ("Executing...")
			
			fftw_execute (fourier_plan);

			data_out [n] = 0.0;
			
			for (int i = 0; i < n; ++i) {
				data_out [i] *= scalar;
			}
			
			TRACE ("Executed.")
		}
		
		inline static std::unique_ptr<plan> make_unique (int i_n, double *i_data_in, double *i_data_out, int *i_flags_ptr = NULL) {
			return std::unique_ptr<plan> (new fftw_cosine (i_n, i_data_in, i_data_out, i_flags_ptr));
		}
	
	private:		
		double scalar;
		fftw_plan fourier_plan;
	};
} /* fft */

#endif /* end of include guard: FFT_HPP_P3TP70YE */
