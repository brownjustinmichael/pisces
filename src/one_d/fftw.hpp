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
		fftw_cosine (int i_n, double *i_data_in, double *i_data_out, int *i_flags_ptr = NULL) : bases::transform (i_n, i_data_in, i_data_out, i_flags_ptr) {
			TRACE ("Instantiating...");
			
			INFO ("FFTW_ESTIMATE = " << FFTW_ESTIMATE);
			
			scalar = 1.0 / sqrt (2.0 * (i_n - 1));
			
			fourier_plan = fftw_plan_r2r_1d (i_n + 1, i_data_in, i_data_out, FFTW_REDFT00, *flags_ptr);
			
			TRACE ("Instantiated.")
		}	
		
		virtual ~fftw_cosine () {}
		
		/*!*******************************************************************
		 * \copydoc bases::transform::execute ()
		 *********************************************************************/
		void execute () {
			TRACE ("Executing...")
			
			fftw_execute (fourier_plan);

			data_out [n] = 0.0;
			
			for (int i = 0; i < n; ++i) {
				data_out [i] *= scalar;
			}
			
			TRACE ("Executed.")
		}
		
		/*!*******************************************************************
		 * \brief Make a unique pointer to a new fftw_cosine object
		 * \copydetails fftw_cosine ()
		 *********************************************************************/
		inline static std::unique_ptr<plan> make_unique (int i_n, double *i_data_in, double *i_data_out, int *i_flags_ptr = NULL) {
			return std::unique_ptr<plan> (new fftw_cosine (i_n, i_data_in, i_data_out, i_flags_ptr));
		}
	
	private:		
		double scalar; //!< The scalar used after the transform (1 / sqrt (2 * (n - 1)))
		fftw_plan fourier_plan; //!< The fftw_plan object to be executed
	};
} /* one_d */

#endif /* end of include guard: FFTW_HPP_P3TP70YE */
