/*!**********************************************************************
 * \file fftw_one_d_cuda.cu
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-16.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <math.h>
#include <fftw3.h>
#include "fftw_one_d_cuda.hpp"
#include "../utils/utils.hpp"
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define PI 3.141592653589793

namespace one_d
{
	namespace cuda
	{
		fftw_cosine::fftw_cosine (bases::element* i_element_ptr, int i_n, int i_name_in, int i_name_out) : 
		bases::transform (i_element_ptr, i_n, i_name_in, i_name_out) {
			scalar = 1.0 / sqrt (2.0 * (n - 1));
			fourier_plan = fftw_plan_r2r_1d (n, data_in, data_out, FFTW_REDFT00, FFTW_ESTIMATE);
		}
		
		void fftw_cosine::execute () {
			bases::transform::execute ();
			if (*flags_ptr & transformed) {
				*flags_ptr &= ~transformed;
			} else {
				*flags_ptr |= transformed;
			}
			
			fftw_execute (fourier_plan);	

			for (int i = 0; i < n + 1; ++i) {
				data_out [i] *= scalar;
			}
		}
	} /* cuda */
} /* one_d */
