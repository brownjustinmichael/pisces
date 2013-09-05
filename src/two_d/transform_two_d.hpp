/*!**********************************************************************
 * \file transform_two_d.hpp
 * /Users/justinbrown/Dropbox/spectral_element/src
 * 
 * Created by Justin Brown on 2013-08-29.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef TRANSFORM_TWO_D_HPP_RAXYBFTC
#define TRANSFORM_TWO_D_HPP_RAXYBFTC

#include "../config.hpp"

namespace bases
{
	template <class datatype>
	class explicit_plan;
} /* bases */

namespace two_d
{
	template <class datatype>
	class explicit_plan : public bases::explicit_plan
	{
	public:
		explicit_plan (int i_n, int i_m, datatype* i_data_in, datatype* i_data_out = NULL) :
		bases::explicit_plan <datatype> (i_n, i_data_in, i_data_out),
		m (i_m) {}
		
		virtual ~explicit_plan () {}
		
		virtual void execute () = 0;
	private:
		int m;
	};
	
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
				
				virtual void execute () {
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
