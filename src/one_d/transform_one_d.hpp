/*!***********************************************************************
 * \file transform_one_d.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-15.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FFTW_HPP_P3TP70YE
#define FFTW_HPP_P3TP70YE

#include <fftw3.h>
#include "../config.hpp"
#include "plan_one_d.hpp"

namespace one_d
{		
	/*!*******************************************************************
	 * An implementation of the transform class using FFTW3.
	 * 
	 * \brief \copybrief bases::transform
	 *********************************************************************/
	template <class datatype>
	class fftw_cosine : public explicit_plan <datatype>
	{
	public:
		/*!*******************************************************************
		 * \copydoc bases::transform::transform ()
		 *********************************************************************/
		fftw_cosine (bases::grid <datatype> &i_grid, datatype* i_data_in, datatype* i_data_out = NULL);
		
		virtual ~fftw_cosine () {}
		
		/*!*******************************************************************
		 * \copydoc bases::transform::execute ()
		 *********************************************************************/
		void execute ();
	
	private:		
		using explicit_plan <datatype>::n;
		using explicit_plan <datatype>::data_in;
		using explicit_plan <datatype>::data_out;
		
		datatype scalar; //!< The scalar used after the transform (1 / sqrt (2 * (n - 1)))
		fftw_plan fourier_plan; //!< The fftw_plan object to be executed
		fftwf_plan fourier_plan_float; //!< The fftw_plan object to be executed
	};
} /* one_d */

#endif /* end of include guard: FFTW_HPP_P3TP70YE */
