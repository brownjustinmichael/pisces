/*!**********************************************************************
 * \file transform.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-15.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cmath>
#include <fftw3.h>
#include "transform.hpp"
#ifdef _MP
#include <omp.h>
#endif
#include "transformer.hpp"

namespace plans
{
	namespace transforms
	{
		template <>
		void horizontal <float>::init (int i_n, int i_m, float* i_data_in, float* i_data_out, int i_flags, int i_threads) {
			n = i_n;
			m = i_m;
			data_in = i_data_in;
			data_out = (i_data_out ? i_data_out : i_data_in);
			flags = i_flags;
		#ifdef _MP
			threads = (i_threads ? i_threads : omp_get_max_threads ());
		#else
			threads = (i_threads ? i_threads : 1);
		#endif
			TRACE ("Initializing...");
			scalar = 1.0 / std::sqrt (n);
			
			// Set up the iodim objects for fftw
			iodim.n = n;
			major_iodim.is = 1;
			major_iodim.os = 1;
			
			plans_float.resize (threads);
			
			if (!(flags & inverse)) {
				// Set up the inverse plan
				iodim.is = m;
				iodim.os = 2 * m;
				int index = 0;
				for (int i = 0; i < threads; ++i) {
					major_iodim.n = m / threads + (i < (m % threads)? 1: 0);
					plans_float [i] = fftwf_plan_guru_split_dft_r2c (1, &iodim, 1, &major_iodim, data_in + index, data_out + index, data_out + m + index, FFTW_MEASURE);
					index += major_iodim.n;
				}
			} else {
				// Set up the forward plan
				iodim.is = 2 * m;
				iodim.os = m;
				int index = 0;
				for (int i = 0; i < threads; ++i) {
					major_iodim.n = m / threads + (i < (m % threads)? 1: 0);
					plans_float [i] = fftwf_plan_guru_split_dft_c2r (1, &iodim, 1, &major_iodim, data_in + index, data_in + m + index, data_out + index, FFTW_MEASURE);
					index += major_iodim.n;
				}
			}
		}
	
		template <>
		horizontal <float>::~horizontal () {
			for (int i = 0; i < threads; ++i) {
				fftwf_destroy_plan (plans_float [i]);
			}
		}
	
		template <>
		void horizontal <double>::init (int i_n, int i_m, double* i_data_in, double* i_data_out, int i_flags, int i_threads) {
			n = i_n;
			m = i_m;
			data_in = i_data_in;
			data_out = (i_data_out ? i_data_out : i_data_in);
			flags = i_flags;
		#ifdef _MP
			threads = (i_threads ? i_threads : omp_get_max_threads ());
		#else
			threads = (i_threads ? i_threads : 1);
		#endif
			TRACE ("Initializing...");
			scalar = 1.0 / std::sqrt (n);
			
			// Set up the iodim objects for fftw
			iodim.n = n;
			major_iodim.is = 1;
			major_iodim.os = 1;
		
			plans.resize (threads);
			
			if (!(flags & inverse)) {
				// Set up the inverse plans
				iodim.is = m;
				iodim.os = 2 * m;
				int index = 0;
				for (int i = 0; i < threads; ++i) {
					major_iodim.n = m / threads + (i < (m % threads)? 1: 0);
					plans [i] = fftw_plan_guru_split_dft_r2c (1, &iodim, 1, &major_iodim, data_in + index, data_out + index, data_out + m + index, FFTW_MEASURE);
					index += major_iodim.n;
				}
			} else {
				// Set up the forward plans
				iodim.is = 2 * m;
				iodim.os = m;
				int index = 0;
				for (int i = 0; i < threads; ++i) {
					major_iodim.n = m / threads + (i < (m % threads)? 1: 0);
					plans [i] = fftw_plan_guru_split_dft_c2r (1, &iodim, 1, &major_iodim, data_in + index, data_in + m + index, data_out + index, FFTW_MEASURE);
					index += major_iodim.n;
				}
			}
		}
		
		template <>
		horizontal <double>::~horizontal () {
			for (int i = 0; i < threads; ++i) {
				fftw_destroy_plan (plans [i]);
			}
		}
	
		template <>
		void horizontal <float>::execute () {
			TRACE ("Executing...");
	// #pragma omp parallel for num_threads (threads)
			for (int i = 0; i < threads; ++i) {
				fftwf_execute (plans_float [i]);
			}
			
			// Switch the horizontal transform flag
			if (*component_flags & transformed_horizontal) {
				*component_flags &= ~transformed_horizontal;
			} else {
				*component_flags |= transformed_horizontal;
			}
		
			linalg::scale (2 * (n / 2 + 1) * m, scalar, data_out);
		}
	
		template <>
		void horizontal <double>::execute () {
			TRACE ("Executing...");
		
	// #pragma omp parallel for num_threads (threads)
			for (int i = 0; i < threads; ++i) {
				fftw_execute (plans [i]);
			}		
			
			// Switch the horizontal transform flag
			if (*component_flags & transformed_horizontal) {
				*component_flags &= ~transformed_horizontal;
			} else {
				*component_flags |= transformed_horizontal;
			}
			
			linalg::scale (2 * (n / 2 + 1) * m, scalar, data_out);
			
			TRACE ("Execution Complete.");
		}
	
		template class horizontal <float>;
		template class horizontal <double>;
	
		template <>
		void vertical <float>::init (int i_n, int i_m, float* i_data_in, float* i_data_out, int i_flags, int i_threads) {
			n = i_n;
			m = i_m;
			data_in = i_data_in;
			data_out = (i_data_out ? i_data_out : i_data_in);
			flags = i_flags;
		#ifdef _MP
			threads = (i_threads ? i_threads : omp_get_max_threads ());
		#else
			threads = (i_threads ? i_threads : 1);
		#endif
			TRACE ("Initializing...");
			scalar = 1.0;
			if (m > 1 && !(flags & ignore_m)) {
				scalar /= std::sqrt (2.0 * (m - 1));
			}
		
			fftwf_r2r_kind kind = FFTW_REDFT00;

			plans_float.resize (threads);
		
			int index = 0;
			for (int i = 0; i < threads; ++i) {
				int nn = (2 * (n / 2 + 1)) / threads + (i < (2 * (n / 2 + 1) % threads)? 1: 0);
				plans_float [i] = fftwf_plan_many_r2r (1, &m, nn, data_in + index * m, NULL, 1, m, data_out + index * m, NULL, 1, m, &kind, FFTW_MEASURE);
				index += nn;
			}
		}
		
		template <>
		vertical <float>::~vertical () {
			fftwf_destroy_plan (plans_float [0]);
		}
		
		template <>
		void vertical <double>::init (int i_n, int i_m, double* i_data_in, double* i_data_out, int i_flags, int i_threads) {
			n = i_n;
			m = i_m;
			data_in = i_data_in;
			data_out = (i_data_out ? i_data_out : i_data_in);
			flags = i_flags;
		#ifdef _MP
			threads = (i_threads ? i_threads : omp_get_max_threads ());
		#else
			threads = (i_threads ? i_threads : 1);
		#endif
			TRACE ("Initializing...");
			scalar = 1.0;
			if (m > 1 && !(flags & ignore_m)) {
				scalar /= std::sqrt (2.0 * (m - 1));
			}
		
			fftw_r2r_kind kind = FFTW_REDFT00;

			plans.resize (threads);
		
			int index = 0;
			for (int i = 0; i < threads; ++i) {
				int nn = (2 * (n / 2 + 1)) / threads + (i < (2 * (n / 2 + 1) % threads)? 1: 0);
				plans [i] = fftw_plan_many_r2r (1, &m, nn, data_in + index * m, NULL, 1, m, data_out + index * m, NULL, 1, m, &kind, FFTW_MEASURE);
				index += nn;
			}
		}
	
		template <>
		vertical <double>::~vertical () {
			fftw_destroy_plan (plans [0]);
		}
	
	
		template <>
		void vertical <float>::execute () {
			TRACE ("Executing...");

			if (m > 1 && !(flags & ignore_m)) {
	// #pragma omp parallel for num_threads (threads)
				for (int i = 0; i < threads; ++i) {
					fftwf_execute (plans_float [i]);
				}
			}
			
			// Switch the vertical transform flag
			if (*component_flags & transformed_vertical) {
				*component_flags &= ~transformed_vertical;
			} else {
				*component_flags |= transformed_vertical;
			}
						
			linalg::scale (2 * (n / 2 + 1) * m, scalar, data_out);
		}
	
		template <>
		void vertical <double>::execute () {
			TRACE ("Executing...");
			
			if (m > 1 && !(flags & ignore_m)) {
	// #pragma omp parallel for num_threads (threads)
				for (int i = 0; i < threads; ++i) {
					fftw_execute (plans [i]);
				}
			}
			
			// Switch the vertical transform flag
			if (*component_flags & transformed_vertical) {
				*component_flags &= ~transformed_vertical;
			} else {
				*component_flags |= transformed_vertical;
			}
		
			linalg::scale (2 * (n / 2 + 1) * m, scalar, data_out);
		
			TRACE ("Execution Complete.");
		}
	
		template class vertical <float>;
		template class vertical <double>;
	} /* transforms */
} /* plans */
