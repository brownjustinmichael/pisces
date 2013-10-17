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
			transform <float>::transform (bases::grid <float> &i_grid_n, bases::grid <float> &i_grid_m, float* i_data_in, float* i_data_out) :
			explicit_plan <float> (i_grid_n, i_grid_m, i_data_in, i_data_out) {
				scalar = 1.0 / std::sqrt (n);
				if (m > 1) {
					scalar /= std::sqrt (2.0 * (m - 1));
				}
				
				iodim_float.n = n;
				iodim_float.is = m;
				iodim_float.os = 2 * m;
				
				major_iodim_float.n = m - 1;
				major_iodim_float.is = 1;
				major_iodim_float.os = 1;
				
				fftwf_r2r_kind kind = FFTW_REDFT00;
				
				x_plan_float = fftwf_plan_guru_split_dft_r2c (1, &iodim, 1, &major_iodim, data_in, data_out, data_out + m, FFTW_ESTIMATE);
				z_plan_float = fftwf_plan_many_r2r (1, &m, n, data_in, NULL, 1, m, data_out, NULL, 1, m, &kind, FFTW_ESTIMATE);
			}
			
			template <>
			transform <double>::transform (bases::grid <double> &i_grid_n, bases::grid <double> &i_grid_m, double* i_data_in, double* i_data_out) :
			explicit_plan <double> (i_grid_n, i_grid_m, i_data_in, i_data_out) {
				scalar = 1.0 / std::sqrt (n);
				if (m > 1) {
					scalar /= std::sqrt (2.0 * (m - 1));
				}
								
				iodim.n = n;
				iodim.is = m;
				iodim.os = 2 * m;
				
				major_iodim.n = m;
				major_iodim.is = 1;
				major_iodim.os = 1;

				fftw_r2r_kind kind = FFTW_REDFT00;
				
				x_plan = fftw_plan_guru_split_dft_r2c (1, &iodim, 1, &major_iodim, data_in, data_out, data_out + m, FFTW_ESTIMATE);
				z_plan = fftw_plan_many_r2r (1, &m, n, data_in, NULL, 1, m, data_out, NULL, 1, m, &kind, FFTW_ESTIMATE);
			}
			
			template <>
			void transform <float>::execute (int &element_flags) {
				TRACE ("Executing...");
	
				// Set up transform
	
				if (m > 1) {
					fftwf_execute (z_plan_float);
				}
				if (n > 1) {
					fftwf_execute (x_plan_float);
				}
				
				// Extract information from transform
	
				for (int i = 0; i < n * m; ++i) {
					data_out [i] *= scalar;
				}
			}
			
			template <>
			void transform <double>::execute (int &element_flags) {
				TRACE ("Executing...");
	
				// Set up transform
	
				if (m > 1) {
					fftw_execute (z_plan);
				}
				if (n > 1) {
					fftw_execute (x_plan);
				}
				
				TRACE ("Scaling...");
				
				// Extract information from transform
	
				for (int i = 0; i < n * m; ++i) {
					data_out [i] *= scalar;
				}
			}
			
			template class transform <float>;
			template class transform <double>;
			
			template <>
			invert <float>::invert (bases::grid <float> &i_grid_n, bases::grid <float> &i_grid_m, float* i_data_in, float* i_data_out) :
			explicit_plan <float> (i_grid_n, i_grid_m, i_data_in, i_data_out) {
				scalar = 1.0 / std::sqrt (n);
				if (m > 1) {
					scalar /= std::sqrt (2.0 * (m - 1));
				}
								
				iodim_float.n = n;
				iodim_float.is = 2 * m;
				iodim_float.os = m;
				
				major_iodim_float.n = m;
				major_iodim_float.is = 1;
				major_iodim_float.os = 1;
				
				fftwf_r2r_kind kind = FFTW_REDFT00;
				
				x_plan_float = fftwf_plan_guru_split_dft_c2r (1, &iodim, 1, &major_iodim, data_in, data_out, data_out + m, FFTW_ESTIMATE);
				z_plan_float = fftwf_plan_many_r2r (1, &m, n, data_in, NULL, 1, m, data_out, NULL, 1, m, &kind, FFTW_ESTIMATE);
			}
			
			template <>
			invert <double>::invert (bases::grid <double> &i_grid_n, bases::grid <double> &i_grid_m, double* i_data_in, double* i_data_out) :
			explicit_plan <double> (i_grid_n, i_grid_m, i_data_in, i_data_out) {
				scalar = 1.0 / std::sqrt (n);
				if (m > 1) {
					scalar /= std::sqrt (2.0 * (m - 1));
				}
								
				iodim.n = n;
				iodim.is = 2 * m;
				iodim.os = m;
				
				major_iodim.n = m;
				major_iodim.is = 1;
				major_iodim.os = 1;

				fftw_r2r_kind kind = FFTW_REDFT00;
				
				x_plan = fftw_plan_guru_split_dft_c2r (1, &iodim, 1, &major_iodim, data_in, data_in + m, data_out, FFTW_ESTIMATE);
				z_plan = fftw_plan_many_r2r (1, &m, n, data_in, NULL, 1, m, data_out, NULL, 1, m, &kind, FFTW_ESTIMATE);
			}
			
			template <>
			void invert <float>::execute (int &element_flags) {
				TRACE ("Executing...");
	
				// Set up invert
	
				if (m > 1) {
					fftwf_execute (z_plan_float);
				}
				if (n > 1) {
					fftwf_execute (x_plan_float);
				}
				
				// Extract information from invert
	
				for (int i = 0; i < n * m; ++i) {
					data_out [i] *= scalar;
				}
			}
			
			template <>
			void invert <double>::execute (int &element_flags) {
				TRACE ("Executing...");
	
				// Set up invert
	
				if (m > 1) {
					fftw_execute (z_plan);
				}
				if (n > 1) {
					fftw_execute (x_plan);
				}
				
				// Extract information from invert
	
				for (int i = 0; i < n * m; ++i) {
					data_out [i] *= scalar;
				}
			}
			
			template class invert <float>;
			template class invert <double>;
		} /* chebyshev */
	} /* fourier */
} /* two_d */