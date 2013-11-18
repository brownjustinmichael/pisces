/*!**********************************************************************
 * \file transform_two_d.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-15.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cmath>
#include <fftw3.h>
#include "transform_two_d.hpp"

namespace two_d
{
	namespace fourier
	{
		namespace chebyshev
		{
			template <>
			horizontal_transform <float>::horizontal_transform (bases::grid <float> &i_grid_n, bases::grid <float> &i_grid_m, float* i_data_in, float* i_data_out, int i_flags) :
			explicit_plan <float> (i_grid_n, i_grid_m, i_data_in, i_data_out),
			flags (i_flags) {
				scalar = 1.0 / std::sqrt (n);
				
				iodim.n = n;
				major_iodim.n = m;
				major_iodim.is = 1;
				major_iodim.os = 1;

				if (!(flags & inverse)) {
					iodim.is = m;
					iodim.os = 2 * m;
				
					x_plan_float = fftwf_plan_guru_split_dft_r2c (1, &iodim, 1, &major_iodim, data_in, data_out, data_out + m, FFTW_ESTIMATE);
				} else {
					iodim.is = 2 * m;
					iodim.os = m;

					x_plan_float = fftwf_plan_guru_split_dft_c2r (1, &iodim, 1, &major_iodim, data_in, data_in + m, data_out, FFTW_ESTIMATE);
				}
			}
			
			template <>
			horizontal_transform <double>::horizontal_transform (bases::grid <double> &i_grid_n, bases::grid <double> &i_grid_m, double* i_data_in, double* i_data_out, int i_flags) :
			explicit_plan <double> (i_grid_n, i_grid_m, i_data_in, i_data_out),
			flags (i_flags) {
				scalar = 1.0 / std::sqrt (n);
								
				iodim.n = n;
				major_iodim.n = m;
				major_iodim.is = 1;
				major_iodim.os = 1;

				if (!(flags & inverse)) {
					iodim.is = m;
					iodim.os = 2 * m;
				
					x_plan = fftw_plan_guru_split_dft_r2c (1, &iodim, 1, &major_iodim, data_in, data_out, data_out + m, FFTW_ESTIMATE);
				} else {
					iodim.is = 2 * m;
					iodim.os = m;

					x_plan = fftw_plan_guru_split_dft_c2r (1, &iodim, 1, &major_iodim, data_in, data_in + m, data_out, FFTW_ESTIMATE);
				}
			}
			
			template <>
			void horizontal_transform <float>::execute (int &element_flags) {
				TRACE ("Executing...");
		
				fftwf_execute (x_plan_float);
					
				for (int i = 0; i < 2 * (n / 2 + 1) * m; ++i) {
					data_out [i] *= scalar;
				}
			}
			
			template <>
			void horizontal_transform <double>::execute (int &element_flags) {
				TRACE ("Executing...");
		
				fftw_execute (x_plan);
				
				for (int i = 0; i < 2 * (n / 2 + 1) * m; ++i) {
					data_out [i] *= scalar;
				}
				
				TRACE ("Execution Complete.");
			}
			
			template class horizontal_transform <float>;
			template class horizontal_transform <double>;
			
			template <>
			vertical_transform <float>::vertical_transform (bases::grid <float> &i_grid_n, bases::grid <float> &i_grid_m, float* i_data_in, float* i_data_out, int i_flags) :
			explicit_plan <float> (i_grid_n, i_grid_m, i_data_in, i_data_out),
			flags (i_flags) {
				scalar = 1.0;
				if (m > 1 && !(flags & ignore_m)) {
					scalar /= std::sqrt (2.0 * (m - 1));
				}
				
				fftwf_r2r_kind kind = FFTW_REDFT00;

				z_plan_float = fftwf_plan_many_r2r (1, &m, n, data_in, NULL, 1, m, data_out, NULL, 1, m, &kind, FFTW_ESTIMATE);
			}
			
			template <>
			vertical_transform <double>::vertical_transform (bases::grid <double> &i_grid_n, bases::grid <double> &i_grid_m, double* i_data_in, double* i_data_out, int i_flags) :
			explicit_plan <double> (i_grid_n, i_grid_m, i_data_in, i_data_out),
			flags (i_flags) {
				scalar = 1.0;
				if (m > 1 && !(flags & ignore_m)) {
					scalar /= std::sqrt (2.0 * (m - 1));
				}

				fftwf_r2r_kind kind = FFTW_REDFT00;

				z_plan = fftw_plan_many_r2r (1, &m, n, data_in, NULL, 1, m, data_out, NULL, 1, m, &kind, FFTW_ESTIMATE);
			}
			
			template <>
			void vertical_transform <float>::execute (int &element_flags) {
				TRACE ("Executing...");
		
				if (m > 1 && !(flags & ignore_m)) {
					fftwf_execute (z_plan_float);
				}
					
				for (int i = 0; i < 2 * (n / 2 + 1) * m; ++i) {
					data_out [i] *= scalar;
				}
			}
			
			template <>
			void vertical_transform <double>::execute (int &element_flags) {
				TRACE ("Executing...");
		
				if (m > 1 && !(flags & ignore_m)) {
					fftw_execute (z_plan);
				}
				
				for (int i = 0; i < 2 * (n / 2 + 1) * m; ++i) {
					data_out [i] *= scalar;
				}
				
				TRACE ("Execution Complete.");
			}
			
			template class vertical_transform <float>;
			template class vertical_transform <double>;
		} /* chebyshev */
	} /* fourier */
} /* two_d */