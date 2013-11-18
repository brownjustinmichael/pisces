/*!**********************************************************************
 * \file solver_two_d.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-14.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cmath>
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
			alpha_0 (grid_m.alpha_0), 
			alpha_n (grid_m.alpha_n),
			n_iterations (i_n_iterations), 
			excess_0 (grid_m.excess_0),
			excess_n (grid_m.excess_n),
			explicit_rhs (i_explicit_rhs), 
			implicit_rhs (i_implicit_rhs) {
				horizontal_plus_matrix.resize (n);
				horizontal_minus_matrix.resize (n);
				factorized_matrix.resize (m * m);
				ipiv.resize (m);
				data_temp.resize (n * m, 0.0);
			}
			
			template <class datatype>
			void solver <datatype>::_factorize () {
				TRACE ("Factorizing...");

				int info;

				for (int i = 0; i < n; ++i) {
					horizontal_plus_matrix [i] = 1.0 + timestep / 2.0 * grid_n.matrix_ptr () [i];
					horizontal_minus_matrix [i] = 1.0 - timestep / 2.0 * grid_n.matrix_ptr () [i];
				}
				
				utils::copy (m * m, grid_m.get_data (0), &factorized_matrix [0]);
		
				utils::add_scaled (m, alpha_0 * timestep / 2.0, grid_m.matrix_ptr () + excess_0, &factorized_matrix [excess_0], m, m);	
				for (int j = excess_0 + 1; j < m - excess_n - 1; ++j) {
					utils::add_scaled (m, timestep / 2.0, grid_m.matrix_ptr () + j, &factorized_matrix [j], m, m);	
				}
				utils::add_scaled (m, alpha_n * timestep / 2.0, grid_m.matrix_ptr () + m - 1 - excess_n, &factorized_matrix [m - 1 - excess_n], m, m);

				utils::matrix_factorize (m, m, &factorized_matrix [0], &ipiv [0], &info);
		
				if (info != 0) {
					ERROR ("Unable to invert matrix");
					throw 0; // For now, kill the program. 
					/*
						TODO Replace this with a more useful exception that can be handled
					*/
				}
			}
			
			template <class datatype>
			void solver <datatype>::execute (int &element_flags) {
				TRACE ("Executing...");

				utils::scale (n * m, 0.0, &data_temp [0]);

				for (int i = 0; i < n; ++i) {
					data_temp [i * m + excess_0] += alpha_0 * timestep / 2.0 * (explicit_rhs [i * m + excess_0] + implicit_rhs [i * m + excess_0]);
					data_temp [i * m + m - 1 - excess_n] += alpha_n * timestep / 2.0 * (explicit_rhs [i * m + m - 1 - excess_n] + implicit_rhs [i * m + m - 1 - excess_n]);
					utils::add_scaled (m - 2 - excess_0 - excess_n, timestep / 2.0, explicit_rhs + i * m + 1 + excess_0, &data_temp [i * m + excess_0 + 1]);
					utils::add_scaled (m - 2 - excess_0 - excess_n, timestep / 2.0, implicit_rhs + i * m + 1 + excess_0, &data_temp [excess_0 + i * m + 1]);
				}
				
				if (element_flags & x_solve) {
					TRACE ("Solving in n direction...");
					for (int j = 0; j < m; ++j) {
						utils::diagonal_multiply (n, 1.0, &horizontal_minus_matrix [0], data_in + j, timestep / 2.0, &data_temp [j], 1, m, m);
						utils::diagonal_solve (n, &horizontal_plus_matrix [0], &data_temp [0] + j, 1, m);
					}
					utils::copy (n * m, &data_temp [0], data_out);
					
					element_flags &= ~x_solve;
					element_flags |= z_solve;
					
				} else if (element_flags & z_solve) {
					TRACE ("Solving in m direction...");
					
					int info;
					utils::add_scaled (n * m, 1.0, data_in, &data_temp [0]);
					
					TRACE ("Beginning matrix solve...");
			
					utils::matrix_solve (m, &factorized_matrix [0], &ipiv [0], &data_temp [0], &info, n, m, m);
		
					TRACE ("Matrix solve complete.");
		
					for (int i = 0; i < n; ++i) {
						for (int j = 0; j < m; ++j) {
							if (std::isnan (data_temp [i * m + j])) {
								ERROR ("Nan detected.");
								info = -1;
							}
						}
					}
		
					if (info != 0) {
						ERROR ("Unable to solve factorized matrix equation");
						throw 0; // For now, kill the program. 
						/*
							TODO Replace this with a more useful exception that can be handled
						*/
					}
			
					TRACE ("Updating...");
					utils::copy (n * m, &data_temp [0], data_out);
					flags |= first_run;
					element_flags |= transformed_vertical;
					element_flags &= ~z_solve;
					element_flags |= x_solve;
		
					TRACE ("Solve complete.")
				}
				TRACE ("Execution complete.");
			}
			
			template class solver <float>;
			template class solver <double>;
		} /* chebyshev */
	} /* fourier */
} /* two_d */