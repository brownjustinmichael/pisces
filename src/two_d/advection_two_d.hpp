/*!**********************************************************************
 * \file advection_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces/src/two_d
 * 
 * Created by Justin Brown on 2013-12-16.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ADVECTION_TWO_D_HPP_GGR0NN1Q
#define ADVECTION_TWO_D_HPP_GGR0NN1Q

#include "plan_two_d.hpp"

namespace two_d
{
	namespace fourier
	{
		// namespace chebyshev
		// {
			template <class datatype>
			class advection : public two_d::real_plan <datatype>
			{
			public:
				advection (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype i_coeff, datatype* i_vel_n, datatype *i_vel_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) :
				real_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out, i_element_flags, i_component_flags),
				coeff (-i_coeff),
				vel_n (i_vel_n),
				vel_m (i_vel_m),
				pos_n (&(grid_n [0])),
				pos_m (&(grid_m [0])) {}
				
				advection (bases::solver <datatype> &i_solver, datatype i_coeff, datatype* i_vel_n, datatype *i_vel_m) :
				real_plan <datatype> (i_solver),
				coeff (-i_coeff),
				vel_n (i_vel_n),
				vel_m (i_vel_m),
				pos_n (&(grid_n [0])),
				pos_m (&(grid_m [0])) {}
				
				virtual ~advection () {}
				
				virtual void execute () {
					#pragma omp parallel for
					for (int j = 1; j < m - 1; ++j) {
						data_out [j] += coeff * (vel_n [j] * (data_in [1 * m + j] - data_in [j]) / (pos_n [1] - pos_n [0]) + vel_m [j] * (data_in [j + 1] - data_in [j - 1]) / (pos_m [j + 1] - pos_m [j - 1]));
						for (int i = 1; i < n - 1; ++i) {
							data_out [i * m + j] += coeff * (vel_n [i * m + j] * (data_in [(i + 1) * m + j] - data_in [(i - 1) * m + j]) / (pos_n [i + 1] - pos_n [i - 1]) + vel_m [i * m + j] * (data_in [i * m + j + 1] - data_in [i * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]));
						}
						data_out [(n - 1) * m + j] += coeff * (vel_n [(n - 1) * m + j] * (data_in [(n - 1) * m + j] - data_in [(n - 2) * m + j]) / (pos_n [n - 1] - pos_n [n - 2]) + vel_m [(n - 1) * m + j] * (data_in [(n - 1) * m + j + 1] - data_in [(n - 1) * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]));
					}
				}
			
			private:
				datatype coeff;
				datatype *vel_n, *vel_m;
				datatype *pos_n, *pos_m;
				using real_plan <datatype>::n;
				using real_plan <datatype>::m;
				using real_plan <datatype>::grid_n;
				using real_plan <datatype>::grid_m;
				using real_plan <datatype>::data_in;
				using real_plan <datatype>::data_out;
			};
			
			template <class datatype>
			class stream_advection : public two_d::real_plan <datatype>
			{
			public:
				stream_advection (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype i_coeff, datatype* i_stream, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) :
				real_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out, i_element_flags, i_component_flags),
				coeff (-i_coeff),
				stream (i_stream),
				pos_n (&(grid_n [0])),
				pos_m (&(grid_m [0])) {}
				
				stream_advection (bases::solver <datatype> &i_solver, datatype i_coeff, datatype *i_stream) :
				real_plan <datatype> (i_solver),
				coeff (-i_coeff),
				stream (i_stream),
				pos_n (&(grid_n [0])),
				pos_m (&(grid_m [0])) {}
				
				virtual ~stream_advection () {}
				
				virtual void execute () {
					for (int j = 1; j < m - 1; ++j) {
						data_out [j] += coeff * ((stream [m + j] - stream [j]) * (data_in [j + 1] - data_in [j - 1]) - (stream [j + 1] - stream [j - 1]) * (data_in [m + j] - data_in [j])) / (pos_n [1] - pos_n [0]) / (pos_m [j + 1] - pos_m [j - 1]);
						for (int i = 1; i < n - 1; ++i) {
							data_out [i * m + j] += coeff * ((stream [(i + 1) * m + j] - stream [(i - 1) * m + j]) * (data_in [i * m + j + 1] - data_in [i * m + j - 1]) - (stream [i * m + j + 1] - stream [i * m + j - 1]) * (data_in [(i + 1) * m + j] - data_in [(i - 1) * m + j])) / (pos_m [j + 1] - pos_m [j - 1]) / (pos_n [i + 1] - pos_n [i - 1]);
						}
						data_out [(n - 1) * m + j] += coeff * ((stream [(n - 1) * m + j] - stream [(n - 2) * m + j]) * (data_in [(n - 1) * m + j + 1] - data_in [(n - 1) * m + j - 1]) - (stream [(n - 1) * m + j + 1] - stream [(n - 1) * m + j - 1]) * (data_in [(n - 1) * m + j] - data_in [(n - 2) * m + j])) / (pos_n [n - 1] - pos_n [n - 2]) / (pos_m [j + 1] - pos_m [j - 1]);
					}
				}
			
			private:
				datatype coeff;
				datatype *stream;
				datatype *pos_n, *pos_m;
				using real_plan <datatype>::n;
				using real_plan <datatype>::m;
				using real_plan <datatype>::grid_n;
				using real_plan <datatype>::grid_m;
				using real_plan <datatype>::data_in;
				using real_plan <datatype>::data_out;
			};
		// } /* chebyshev */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: ADVECTION_TWO_D_HPP_GGR0NN1Q */