/*!**********************************************************************
 * \file solver_two_d.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-14.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "solver_two_d.hpp"

namespace two_d
{
	namespace fourier
	{
		namespace chebyshev
		{
			template <class datatype>
			solver::_factorize () {
				for (int i = 0; i < n; ++i) {
					horizontal_plus_matrix [i] = 1.0 + timestep * alpha * matrix_n [i];
					horizontal_minus_matrix [i] = 1.0 - timestep * (1.0 - alpha) * matrix_n [i];
				}
			}
			
			template <class datatype>
			solver::execute (int &element_flags) {
				utils::copy (n * m, explicit_rhs, &data_temp [0]);
				utils::add_scaled (n * m, 1.0, implicit_rhs, &data_temp [0]);
				for (int j = 0; j < m; ++j) {
					utils::diagonal_multiply (n, 1.0, &horizontal_minus_matrix [0], data_in + j, timestep, &data_temp [0] + j, 1, m, m);
					utils::diagonal_solve (n, &horizontal_plus_matrix [0], &data_temp [0] + j, 1, m);
				}
				utils::copy (n * m, &data_temp [0], data_out);
			}
			
			template class solver <float>
			template class solver <double>
		} /* chebyshev */
	} /* fourier */
} /* two_d */