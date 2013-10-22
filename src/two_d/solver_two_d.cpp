/*!**********************************************************************
 * \file solver_two_d.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-14.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "../utils/utils.hpp"
#include "../utils/solver_utils.hpp"
#include "solver_two_d.hpp"

namespace two_d
{
	namespace fourier
	{
		namespace chebyshev
		{
			template <class datatype>
			solver <datatype>:: solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, bases::messenger* i_messenger_ptr, int i_n_iterations, datatype& i_timestep, datatype* i_data_in, datatype* i_explicit_rhs, datatype* i_implicit_rhs, datatype* i_data_out, int i_flags) : 
			bases::solver <datatype> (i_flags),
			explicit_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out), 
			messenger_ptr (i_messenger_ptr), 
			timestep (i_timestep), 
			alpha_0 (grid_n.alpha_0), 
			alpha_n (grid_n.alpha_n),
			n_iterations (i_n_iterations), 
			explicit_rhs (i_explicit_rhs), 
			implicit_rhs (i_implicit_rhs) {
				horizontal_plus_matrix.resize (n);
				horizontal_minus_matrix.resize (n);
				data_temp.resize (n * m);
			}
			
			template <class datatype>
			void solver <datatype>::_factorize () {
				for (int i = 0; i < n; ++i) {
					horizontal_plus_matrix [i] = 1.0 + timestep * grid_n.matrix_ptr () [i];
					horizontal_minus_matrix [i] = 1.0 - timestep * grid_n.matrix_ptr () [i];
				}
			}
			
			template <class datatype>
			void solver <datatype>::execute (int &element_flags) {
				TRACE ("Executing...");
				
				utils::copy (n * m, explicit_rhs, &data_temp [0]);
				utils::add_scaled (n * m, 1.0, implicit_rhs, &data_temp [0]);
				for (int j = 0; j < m; ++j) {
					utils::diagonal_multiply (n, 1.0, &horizontal_minus_matrix [0], data_in + j, timestep, &data_temp [0] + j, 1, m, m);
					utils::diagonal_solve (n, &horizontal_plus_matrix [0], &data_temp [0] + j, 1, m);
				}
				utils::copy (n * m, &data_temp [0], data_out);
				TRACE ("Execution complete.");
			}
			
			template class solver <float>;
			template class solver <double>;
		} /* chebyshev */
	} /* fourier */
} /* two_d */