/*!**********************************************************************
 * \file transform_two_d.cpp
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
	template <>
	horizontal_transform <float>::horizontal_transform (int n, int m, float* i_data_in, float* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads) :
	plans::plan <float> (i_element_flags, i_component_flags),
	n (n), 
	m (m), 
	data_in (i_data_in),
	data_out (i_data_out ? i_data_out : i_data_in),
	flags (i_flags),
#ifdef _MP
	threads (i_threads ? i_threads : omp_get_max_threads ()) {
#else
	threads (i_threads ? i_threads : 1) {
#endif
		TRACE ("Initializing...");
		scalar = 1.0 / std::sqrt (n);
						
		iodim.n = n;
		major_iodim.is = 1;
		major_iodim.os = 1;
		
		plans_float.resize (threads);
		
		if (!(flags & inverse)) {
			iodim.is = m;
			iodim.os = 2 * m;
			int index = 0;
			for (int i = 0; i < threads; ++i) {
				major_iodim.n = m / threads + (i < (m % threads)? 1: 0);
				plans_float [i] = fftwf_plan_guru_split_dft_r2c (1, &iodim, 1, &major_iodim, data_in + index, data_out + index, data_out + m + index, FFTW_MEASURE);
				index += major_iodim.n;
			}
		} else {
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
	horizontal_transform <float>::horizontal_transform (plans::grid <float> &i_grid_n, plans::grid <float> &i_grid_m, float* i_data_in, float* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads) :
	plans::plan <float> (i_element_flags, i_component_flags),
	n (i_grid_n.get_n ()), 
	m (i_grid_m.get_n ()), 
	data_in (i_data_in),
	data_out (i_data_out ? i_data_out : i_data_in),
	flags (i_flags),
#ifdef _MP
	threads (i_threads ? i_threads : omp_get_max_threads ()) {
#else
	threads (i_threads ? i_threads : 1) {
#endif
		TRACE ("Initializing...");
		scalar = 1.0 / std::sqrt (n);
						
		iodim.n = n;
		major_iodim.is = 1;
		major_iodim.os = 1;
		
		plans_float.resize (threads);
		
		if (!(flags & inverse)) {
			iodim.is = m;
			iodim.os = 2 * m;
			int index = 0;
			for (int i = 0; i < threads; ++i) {
				major_iodim.n = m / threads + (i < (m % threads)? 1: 0);
				plans_float [i] = fftwf_plan_guru_split_dft_r2c (1, &iodim, 1, &major_iodim, data_in + index, data_out + index, data_out + m + index, FFTW_MEASURE);
				index += major_iodim.n;
			}
		} else {
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
	horizontal_transform <float>::~horizontal_transform () {
		for (int i = 0; i < threads; ++i) {
			fftwf_destroy_plan (plans_float [i]);
		}
	}
	
	template <>
	horizontal_transform <double>::horizontal_transform (int n, int m, double* i_data_in, double* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads) :
	plans::plan <double> (i_element_flags, i_component_flags),
	n (n), 
	m (m), 
	data_in (i_data_in),
	data_out (i_data_out ? i_data_out : i_data_in),
	flags (i_flags),
#ifdef _MP
	threads (i_threads ? i_threads : omp_get_max_threads ()) {
#else
	threads (i_threads ? i_threads : 1) {
#endif
		TRACE ("Initializing...");
		scalar = 1.0 / std::sqrt (n);
						
		iodim.n = n;
		major_iodim.is = 1;
		major_iodim.os = 1;
		
		plans.resize (threads);
		
		if (!(flags & inverse)) {
			iodim.is = m;
			iodim.os = 2 * m;
			int index = 0;
			for (int i = 0; i < threads; ++i) {
				major_iodim.n = m / threads + (i < (m % threads)? 1: 0);
				plans [i] = fftw_plan_guru_split_dft_r2c (1, &iodim, 1, &major_iodim, data_in + index, data_out + index, data_out + m + index, FFTW_MEASURE);
				index += major_iodim.n;
			}
		} else {
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
	horizontal_transform <double>::horizontal_transform (plans::grid <double> &i_grid_n, plans::grid <double> &i_grid_m, double* i_data_in, double* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads) :
	plans::plan <double> (i_element_flags, i_component_flags),
	n (i_grid_n.get_n ()), 
	m (i_grid_m.get_n ()), 
	data_in (i_data_in),
	data_out (i_data_out ? i_data_out : i_data_in),
	flags (i_flags),
#ifdef _MP
	threads (i_threads ? i_threads : omp_get_max_threads ()) {
#else
	threads (i_threads ? i_threads : 1) {
#endif
		TRACE ("Initializing...");
		scalar = 1.0 / std::sqrt (n);
						
		iodim.n = n;
		major_iodim.is = 1;
		major_iodim.os = 1;
		
		plans.resize (threads);
		
		if (!(flags & inverse)) {
			iodim.is = m;
			iodim.os = 2 * m;
			int index = 0;
			for (int i = 0; i < threads; ++i) {
				major_iodim.n = m / threads + (i < (m % threads)? 1: 0);
				plans [i] = fftw_plan_guru_split_dft_r2c (1, &iodim, 1, &major_iodim, data_in + index, data_out + index, data_out + m + index, FFTW_MEASURE);
				index += major_iodim.n;
			}
		} else {
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
	horizontal_transform <double>::~horizontal_transform () {
		for (int i = 0; i < threads; ++i) {
			fftw_destroy_plan (plans [i]);
		}
	}
	
	template <>
	void horizontal_transform <float>::execute () {
		TRACE ("Executing...");
// #pragma omp parallel for num_threads (threads)
		for (int i = 0; i < threads; ++i) {
			fftwf_execute (plans_float [i]);
		}				
		
		if (*component_flags & transformed_horizontal) {
			*component_flags &= ~transformed_horizontal;
		} else {
			*component_flags |= transformed_horizontal;
			// for (int i = 4 * (n / 2 + 1) / 3; i < 2 * (n / 2 + 1); ++i) {
			// 	for (int j = 0; j < m; ++j) {
			// 		data_out [i * m + j] *= 0.0;
			// 	}
			// }
		}
			
		for (int i = 0; i < 2 * (n / 2 + 1) * m; ++i) {
			data_out [i] *= scalar;
		}
	}
	
	template <>
	void horizontal_transform <double>::execute () {
		TRACE ("Executing...");
		
		// std::stringstream debug;
		// for (int i = 0; i < 2 * (n / 2 + 1); ++i) {
		// 	debug << data_in [i * m + m - 1] << " ";
		// }
		// DEBUG ("IN " << debug.str ());
		// debug.str ("");

// #pragma omp parallel for num_threads (threads)
		for (int i = 0; i < threads; ++i) {
			fftw_execute (plans [i]);
		}		
		
		if (*component_flags & transformed_horizontal) {
			*component_flags &= ~transformed_horizontal;
		} else {
			*component_flags |= transformed_horizontal;
			// for (int i = 4 * (n / 2 + 1) / 3; i < 2 * (n / 2 + 1); ++i) {
			// 	for (int j = 0; j < m; ++j) {
			// 		data_in [i * m + j] = 0.0;
			// 	}
			// }
		}
	
		for (int i = 0; i < 2 * (n / 2 + 1) * m; ++i) {
			data_out [i] *= scalar;
		}
		// for (int i = 0; i < 2 * (n / 2 + 1); ++i) {
		// 	debug << data_out [i * m + m - 1] << " ";
		// }
		// DEBUG ("OUT " << debug.str ());
		TRACE ("Execution Complete.");
	}
	
	template class horizontal_transform <float>;
	template class horizontal_transform <double>;
	
	template <>
	vertical_transform <float>::vertical_transform (int n, int m, float* i_data_in, float* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads) :
	plans::plan <float> (i_element_flags, i_component_flags),
	n (n), 
	m (m), 
	data_in (i_data_in),
	data_out (i_data_out ? i_data_out : i_data_in),
	flags (i_flags),
#ifdef _MP
	threads (i_threads ? i_threads : omp_get_max_threads ()) {
#else
	threads (i_threads ? i_threads : 1) {
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
	vertical_transform <float>::vertical_transform (plans::grid <float> &i_grid_n, plans::grid <float> &i_grid_m, float* i_data_in, float* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads) :
	plans::plan <float> (i_element_flags, i_component_flags),
	n (i_grid_n.get_n ()), 
	m (i_grid_m.get_n ()), 
	data_in (i_data_in),
	data_out (i_data_out ? i_data_out : i_data_in),
	flags (i_flags),
#ifdef _MP
	threads (i_threads ? i_threads : omp_get_max_threads ()) {
#else
	threads (i_threads ? i_threads : 1) {
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
	vertical_transform <float>::~vertical_transform () {
		fftwf_destroy_plan (plans_float [0]);
	}
	
	template <>
	vertical_transform <double>::vertical_transform (int n, int m, double* i_data_in, double* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads) :
	plans::plan <double> (i_element_flags, i_component_flags),
	n (n), 
	m (m), 
	data_in (i_data_in),
	data_out (i_data_out ? i_data_out : i_data_in),
	flags (i_flags),
#ifdef _MP
	threads (i_threads ? i_threads : omp_get_max_threads ()) {
#else
	threads (i_threads ? i_threads : 1) {
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
	vertical_transform <double>::vertical_transform (plans::grid <double> &i_grid_n, plans::grid <double> &i_grid_m, double* i_data_in, double* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads) :
	plans::plan <double> (i_element_flags, i_component_flags),
	n (i_grid_n.get_n ()), 
	m (i_grid_m.get_n ()), 
	data_in (i_data_in),
	data_out (i_data_out ? i_data_out : i_data_in),
	flags (i_flags),
#ifdef _MP
	threads (i_threads ? i_threads : omp_get_max_threads ()) {
#else
	threads (i_threads ? i_threads : 1) {
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
	vertical_transform <double>::~vertical_transform () {
		fftw_destroy_plan (plans [0]);
	}
	
	
	template <>
	void vertical_transform <float>::execute () {
		TRACE ("Executing...");

		if (m > 1 && !(flags & ignore_m)) {
// #pragma omp parallel for num_threads (threads)
			for (int i = 0; i < threads; ++i) {
				fftwf_execute (plans_float [i]);
			}
		}
		
		if (*component_flags & transformed_vertical) {
			*component_flags &= ~transformed_vertical;
		} else {
			*component_flags |= transformed_vertical;
		}
							
		for (int i = 0; i < 2 * (n / 2 + 1) * m; ++i) {
			data_out [i] *= scalar;
		}
	}
	
	template <>
	void vertical_transform <double>::execute () {
		TRACE ("Executing...");
		// DEBUG ("VERTI");
		// std::stringstream debug;
		// for (int j = 0; j < m; ++j) {
		// 	for (int i = 0; i < 2 * (n / 2 + 1); ++i) {
		// 		debug << data_in [i * m + j] << " ";
		// 	}
		// 	DEBUG ("IN " << debug.str ());
		// 	debug.str ("");
		// }

		if (m > 1 && !(flags & ignore_m)) {
// #pragma omp parallel for num_threads (threads)
			for (int i = 0; i < threads; ++i) {
				fftw_execute (plans [i]);
			}
		}

		if (*component_flags & transformed_vertical) {
			*component_flags &= ~transformed_vertical;
			// DEBUG ("A");
		} else {
			*component_flags |= transformed_vertical;
			// DEBUG ("B");
		}
		
		for (int i = 0; i < 2 * (n / 2 + 1) * m; ++i) {
			data_out [i] *= scalar;
		}
		
		// for (int j = 0; j < m; ++j) {
		// 	for (int i = 0; i < 2 * (n / 2 + 1); ++i) {
		// 		debug << data_out [i * m + j] << " ";
		// 	}
		// 	DEBUG ("OUT " << debug.str ());
		// 	debug.str ("");
		// }
		
		TRACE ("Execution Complete.");
	}
	
	template class vertical_transform <float>;
	template class vertical_transform <double>;
} /* plans */