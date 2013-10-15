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


namespace two_d
{
	namespace chebyshev
	{
		namespace fourier
		{
			template <class datatype>
			class transform
			{
			public:
				transform (int i_n, int i_m, datatype* i_data_in, datatype* i_data_out = NULL) :
				explicit_plan (i_n, i_m, i_data_in, i_data_out) {
					scalar = 1.0 / std::sqrt (2.0 * (n - 1)) / std::sqrt (2.0 * (m - 1));
					fourier_plan = fftw_plan_r2c_2d (n, m, data_in, data_out, FFTW_ESTIMATE);
				}
				
				virtual ~transform () {}
				
				virtual void execute (int &element_flags) {
					TRACE ("Executing...");
		
					// Set up transform
		
					fftw_execute (fourier_plan);
					
					// Extract information from transform
		
					for (int i = 0; i < n * m; ++i) {
						data_out [i] *= scalar;
					}
				}
			
			protected:
				using explicit_plan <datatype>::n;
				using explicit_plan <datatype>::m;
				using explicit_plan <datatype>::data_in;
				using explicit_plan <datatype>::data_out;
				
				datatype scalar;
				fftw_plan fourier_plan;
			};
		} /* fourier */
	} /* chebyshev */
} /* two_d */

#endif /* end of include guard: TRANSFORM_TWO_D_HPP_RAXYBFTC */
