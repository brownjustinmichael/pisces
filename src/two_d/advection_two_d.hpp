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
		namespace chebyshev
		{
			template <class datatype>
			class advection : public two_d::explicit_plan <datatype>
			{
			public:
				advection (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype i_coeff, datatype* i_vel_n, datatype *i_vel_m, datatype* i_data_in = NULL, datatype *i_data_out = NULL) :
				explicit_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out),
				coeff (-i_coeff),
				vel_n (i_vel_n),
				vel_m (i_vel_m),
				pos_n (&(i_grid_n.position ())),
				pos_m (&(i_grid_m.position ())) {}
				
				virtual ~advection () {}
				
				virtual void execute (bases::flags &element_flags) {
					for (int j = 1; j < m - 1; ++j) {
						for (int i = 1; i < n - 1; ++i) {
							data_out [i * m + j] += coeff * (vel_n [i * m + j] * (data_in [(i + 1) * m + j] - data_in [(i - 1) * m + j]) / (pos_n [i + 1] - pos_n [i - 1]) + vel_m [i * m + j] * (data_in [i * m + j + 1] - data_in [i * m + j - 1]) / (pos_m [j + 1] - pos_n [j - 1]));
						}
					}
				}
			
			private:
				datatype coeff;
				datatype *vel_n, *vel_m;
				datatype *pos_n, *pos_m;
				using explicit_plan <datatype>::n;
				using explicit_plan <datatype>::m;
				using explicit_plan <datatype>::data_in;
				using explicit_plan <datatype>::data_out;
			};
			
			template <class datatype>
			class stream_advection : public two_d::explicit_plan <datatype>
			{
			public:
				stream_advection (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype i_coeff, datatype *i_stream, datatype* i_data_in = NULL, datatype *i_data_out = NULL) :
				explicit_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out),
				coeff (-i_coeff),
				stream (i_stream),
				pos_n (&(i_grid_n.position ())),
				pos_m (&(i_grid_m.position ())) {}
				
				virtual ~stream_advection () {}
				
				virtual void execute (bases::flags &element_flags) {
					for (int j = 1; j < m - 1; ++j) {
						for (int i = 1; i < n - 1; ++i) {
							data_out [i * m + j] += coeff * ((stream [(i + 1) * m + j] - stream [(i - 1) * m + j]) * (data_in [i * m + j + 1] - data_in [i * m + j - 1]) + (stream [i * m + j + 1] - stream [i * m + j - 1]) * (data_in [(i + 1) * m + j] - data_in [(i - 1) * m + j])) / (pos_m [j + 1] - pos_n [j - 1]) / (pos_n [i + 1] - pos_n [i - 1]);
						}
					}
				}
			
			private:
				datatype coeff;
				datatype *stream;
				datatype *pos_n, *pos_m;
				using explicit_plan <datatype>::n;
				using explicit_plan <datatype>::m;
				using explicit_plan <datatype>::data_in;
				using explicit_plan <datatype>::data_out;
			};
		} /* chebyshev */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: ADVECTION_TWO_D_HPP_GGR0NN1Q */
