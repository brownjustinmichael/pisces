/*!***********************************************************************
 * \file one_d/fftw.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-15.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FFTW_HPP_P3TP70YE
#define FFTW_HPP_P3TP70YE

#include <fftw3.h>
#include "../config.hpp"
#include "../bases/transform.hpp"

namespace one_d
{		
	/*!*******************************************************************
	 * An implementation of the transform class using FFTW3.
	 * 
	 * \brief \copybrief bases::transform
	 *********************************************************************/
	class fftw_cosine : public bases::transform
	{
	public:
		/*!*******************************************************************
		 * \copydoc bases::transform::transform ()
		 *********************************************************************/
		fftw_cosine (int i_n, double *i_data_in, double *i_data_out = NULL, int *i_flags_ptr = NULL, int i_logger = -1) : bases::transform (i_n, i_data_in, i_data_out, i_flags_ptr, i_logger) {
			init (i_n, i_data_in, i_data_out, i_flags_ptr);
		}
		
		fftw_cosine (int i_n, double& i_data_in, double& i_data_out, int *i_flags_ptr = NULL, int i_logger = -1) : bases::transform (i_n, &i_data_in, &i_data_out, i_flags_ptr) {
			init (i_n, &i_data_in, &i_data_out, i_flags_ptr);
		}
		
		fftw_cosine (int i_n, double& i_data_in, int *i_flags_ptr = NULL, int i_logger = -1) : bases::transform (i_n, &i_data_in, NULL, i_flags_ptr, i_logger) {
			init (i_n, &i_data_in);
		}
		
		virtual ~fftw_cosine () {}
		
		/*!*******************************************************************
		 * \copydoc bases::transform::execute ()
		 *********************************************************************/
		void execute () {
			TRACE (logger, "Executing...");
			
			if (*flags_ptr & transformed) {
				*flags_ptr &= ~transformed;
			} else {
				*flags_ptr |= transformed;
			}
			
			fftw_execute (fourier_plan);
			
			if (*flags_ptr & transformed) {
				// data_out [n - 1] = 0.0;
			}
			
			for (int i = 0; i < n; ++i) {
				data_out [i] *= scalar;
			}
			
			TRACE (logger, "Executed.");
		}
	
	private:		
		double scalar; //!< The scalar used after the transform (1 / sqrt (2 * (n - 1)))
		fftw_plan fourier_plan; //!< The fftw_plan object to be executed
		int runs;
		
		inline void init (int i_n, double *i_data_in, double *i_data_out = NULL, int *i_flags_ptr = NULL) {
			TRACE (logger, "Instantiating...");
						
			scalar = 1.0 / sqrt (2.0 * (i_n - 1));
			
			fourier_plan = fftw_plan_r2r_1d (i_n, data_in, data_out, FFTW_REDFT00, *flags_ptr);
			
			TRACE (logger, "Instantiated.")
		}
	};
} /* one_d */

#endif /* end of include guard: FFTW_HPP_P3TP70YE */
