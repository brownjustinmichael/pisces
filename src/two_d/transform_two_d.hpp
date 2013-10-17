/*!**********************************************************************
 * \file transform_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-29.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef TRANSFORM_TWO_D_HPP_RAXYBFTC
#define TRANSFORM_TWO_D_HPP_RAXYBFTC

#include "../config.hpp"
#include "plan_two_d.hpp"
#include <fftw3.h>

namespace two_d
{
	namespace fourier
	{
		namespace chebyshev
		{
			template <class datatype>
			class transform : public explicit_plan <datatype>
			{
			public:
				transform (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out = NULL);
				
				virtual ~transform () {}
				
				virtual void execute (int &element_flags);
			
			protected:
				using explicit_plan <datatype>::n;
				using explicit_plan <datatype>::m;
				using explicit_plan <datatype>::data_in;
				using explicit_plan <datatype>::data_out;
				
				datatype scalar;
				fftw_plan x_plan;
				fftw_plan z_plan;
				fftwf_plan x_plan_float;
				fftwf_plan z_plan_float;
				fftw_iodim major_iodim;
				fftw_iodim iodim;
				fftwf_iodim major_iodim_float;
				fftwf_iodim iodim_float;

			};
			
			template <class datatype>
			class invert : public explicit_plan <datatype>
			{
			public:
				invert (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out = NULL);
				
				virtual ~invert () {}
				
				virtual void execute (int &element_flags);
			
			protected:
				using explicit_plan <datatype>::n;
				using explicit_plan <datatype>::m;
				using explicit_plan <datatype>::data_in;
				using explicit_plan <datatype>::data_out;
				
				datatype scalar;
				fftw_plan x_plan;
				fftw_plan z_plan;
				fftwf_plan x_plan_float;
				fftwf_plan z_plan_float;
				fftw_iodim major_iodim;
				fftw_iodim iodim;
				fftwf_iodim major_iodim_float;
				fftwf_iodim iodim_float;

			};
		} /* chebyshev */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: TRANSFORM_TWO_D_HPP_RAXYBFTC */
