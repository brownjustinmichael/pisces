/*!***********************************************************************
 * \file fftw_one_d.hpp
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

namespace bases
{
	class element;
} /* bases */

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
		fftw_cosine (bases::element* i_element_ptr, int i_n, int i_name_in, int i_name_out = null) : bases::transform (i_element_ptr, i_n, i_name_in, i_name_out) {
			scalar = 1.0 / sqrt (2.0 * (n - 1));
			fourier_plan = fftw_plan_r2r_1d (n, data_in, data_out, FFTW_REDFT00, FFTW_ESTIMATE);
		}
		
		virtual ~fftw_cosine () {}
		
		/*!*******************************************************************
		 * \copydoc bases::transform::execute ()
		 *********************************************************************/
		void execute () {
			TRACE ("Executing...");
			
			bases::transform::execute ();
			
			if (*flags_ptr & transformed) {
				*flags_ptr &= ~transformed;
			} else {
				*flags_ptr |= transformed;
			}
			
			fftw_execute (fourier_plan);
			
			for (int i = 0; i < n; ++i) {
				data_out [i] *= scalar;
			}
			
		}
	
	private:		
		double scalar; //!< The scalar used after the transform (1 / sqrt (2 * (n - 1)))
		fftw_plan fourier_plan; //!< The fftw_plan object to be executed
	};
} /* one_d */

#endif /* end of include guard: FFTW_HPP_P3TP70YE */
