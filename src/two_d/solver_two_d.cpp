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
#include "../utils/interpolate.hpp"
#include "../utils/block_solver.hpp"
#include "solver_two_d.hpp"

namespace two_d
{
	namespace fourier
	{
		namespace chebyshev
		{
			template <class datatype>
			solver <datatype>:: solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, bases::messenger* i_messenger_ptr, datatype& i_timestep, datatype& i_alpha_0, datatype& i_alpha_n, datatype* i_data_in, datatype* i_explicit_rhs, datatype* i_implicit_rhs, datatype* i_data_out, int i_flags) : 
			bases::solver <datatype> (i_flags),
			explicit_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out), 
			messenger_ptr (i_messenger_ptr), 
			timestep (i_timestep), 
			alpha_0 (i_alpha_0), 
			alpha_n (i_alpha_n), 
			positions (&(grid_m.position ())),
			excess_0 (grid_m.excess_0), 
			excess_n (grid_m.excess_n),
			explicit_rhs (i_explicit_rhs),
			implicit_rhs (i_implicit_rhs),
			default_matrix (grid_m.get_data (0)), 
			matrix (grid_m.matrix_ptr ()) {
				horizontal_plus_matrix.resize (n);
				horizontal_minus_matrix.resize (n);
				values_0.resize (n);
				values_n.resize (n);
				if (messenger_ptr->get_id () - 1 >= 0) {
					messenger_ptr->template send <int> (1, &excess_0, messenger_ptr->get_id () - 1, 0);
					messenger_ptr->template recv <int> (1, &ex_excess_0, messenger_ptr->get_id () - 1, 0);
					positions_0.resize (ex_excess_0);
					messenger_ptr->template send <datatype> (excess_0, positions, messenger_ptr->get_id () - 1, 0);
					messenger_ptr->template recv <datatype> (ex_excess_0, &positions_0 [0], messenger_ptr->get_id () - 1, 0);
					ntop = 1;
				} else {
					ex_excess_0 = 0;
					ntop = 0;
				}
				if (messenger_ptr->get_id () + 1 < messenger_ptr->get_np ()) {
					messenger_ptr->template send <int> (1, &excess_n, messenger_ptr->get_id () + 1, 0);
					messenger_ptr->template recv <int> (1, &ex_excess_n, messenger_ptr->get_id () + 1, 0);
					positions_n.resize (ex_excess_n);
					messenger_ptr->template send <datatype> (excess_n, &(positions [m - excess_n]), messenger_ptr->get_id () + 1, 0);
					messenger_ptr->template recv <datatype> (ex_excess_n, &positions_n [0], messenger_ptr->get_id () + 1, 0);
					nbot = 1;
				} else {
					ex_excess_n = 0;
					nbot = 0;
				}
				int ns0 = excess_0 + ex_excess_0 + ntop * 2;
				if (messenger_ptr->get_id () == 0) {
					ns.resize (messenger_ptr->get_np ());
					messenger_ptr->template gather <int> (1, &ns0, &ns [0]);
					int ntot = 0;
					for (int i = 0; i < messenger_ptr->get_np (); ++i) {
						ntot += ns [i];
					}
					boundary_matrix.resize (ntot * ntot);
					bipiv.resize (ntot);
				} else {
					messenger_ptr->template gather <int> (1, &ns0, NULL);
					boundary_matrix.resize ((excess_0 + ex_excess_0 + excess_n + ex_excess_n + 2 * (nbot + ntop)) * (excess_0 + ex_excess_0 + excess_n + ex_excess_n + 2 * (nbot + ntop)));
				}
				factorized_matrix.resize ((m + ex_excess_0 + ex_excess_n + nbot + ntop) * (m + ex_excess_0 + ex_excess_n + nbot + ntop), 0.0);
				ipiv.resize (m); // Should be n - ntop - nbot - excess_0 - excess_n
				data_temp.resize ((m + ex_excess_0 + ex_excess_n + nbot + ntop) * n);
			}
			
			template <class datatype>
			void solver <datatype>::_factorize () {
				int info, lda = m + ex_excess_0 + ex_excess_n + nbot + ntop;
				TRACE ("Factorizing...");
				
				for (int i = 0; i < n; ++i) {
					horizontal_plus_matrix [i] = 1.0 + timestep / 2.0 * grid_n.matrix_ptr () [i];
					horizontal_minus_matrix [i] = 1.0 - timestep / 2.0 * grid_n.matrix_ptr () [i];
				}
						
				utils::matrix_copy (m, m, default_matrix, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1)], m, lda);

				utils::matrix_add_scaled (m - excess_n - excess_0 - 2, m, timestep, matrix + excess_0 + 1, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + 1 + excess_0], m, lda);
				if (ntop != 0) {
					utils::matrix_add_scaled (ntop, m, alpha_0 * timestep, matrix + excess_0, &factorized_matrix [(ntop + ex_excess_0) * lda], m, m + ex_excess_0 + ex_excess_n + ntop + nbot);
					utils::interpolate (ex_excess_0, m, m, timestep, positions, matrix, &positions_0 [0], &factorized_matrix [(ntop + ex_excess_0) * lda + ntop], m, lda);
					utils::matrix_add_scaled (ntop, m, alpha_0 * timestep, matrix + excess_0, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + excess_0], m, m + ex_excess_0 + ex_excess_n + ntop + nbot);
				}
				if (nbot != 0) {
					utils::matrix_add_scaled (nbot, m, alpha_n * timestep, matrix + m - nbot - excess_n, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + m - nbot - excess_n], m, lda);
					utils::interpolate (ex_excess_n, m, m, timestep, positions, matrix, &positions_n [0], &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + m], m, lda);
					utils::matrix_add_scaled (nbot, m, alpha_n * timestep, matrix + m - nbot - excess_n, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + m + ex_excess_n], m, lda);
				}

				utils::p_block_matrix_factorize (messenger_ptr->get_id (), messenger_ptr->get_np (), m - excess_0 - excess_n - ntop - nbot, excess_0 + ex_excess_0 + 2 * ntop, excess_n + ex_excess_n + 2 * nbot, &factorized_matrix [0], &ipiv [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, lda, sqrt ((int) boundary_matrix.size ()));
				
				if (info != 0) {
					ERROR ("Unable to invert matrix");
					throw 0; // For now, kill the program. 
					/*
						TODO Replace this with a more useful exception that can be handled
					*/
				}
				TRACE ("Done.");
			}
			
			template <class datatype>
			void solver <datatype>::execute (int &element_flags) {
				int info, lda = m + ex_excess_0 + ex_excess_n + nbot + ntop;
				TRACE ("Executing solve...");
				utils::scale (n * lda, 0.0, &data_temp [0]);
						
				if (!(flags & first_run)) {
					utils::copy (n, &data_in [0], &values_0 [0], m);
					utils::copy (n, &data_in [m - 1], &values_n [0], m);
					flags |= first_run;
				}
				
										
				utils::matrix_add_scaled (m - excess_0 - excess_n, n, timestep / 2.0, implicit_rhs + excess_0, &data_temp [ex_excess_0 + ntop + excess_0], m, lda);
				utils::matrix_add_scaled (m - excess_0 - excess_n, n, timestep / 2.0, explicit_rhs + excess_0, &data_temp [ex_excess_0 + ntop + excess_0], m, lda);
								
				if (ntop != 0) {
					utils::interpolate (ex_excess_0, n, m - excess_0 - excess_n, 1.0, positions + excess_0, &data_temp [ex_excess_0 + ntop + excess_0], &positions_0 [0], &data_temp [ntop], lda, lda);
					utils::scale (n, alpha_0, &data_temp [ntop + ex_excess_0 + excess_0], lda);
					utils::copy (n, &data_temp [0], &data_temp [ntop + ex_excess_0 + excess_0], lda, lda);
				} else {
					utils::copy (n, &values_0 [0], &data_temp [0], 1, lda);
				}
				if (nbot != 0) {
					utils::interpolate (ex_excess_n, n, m - excess_0 - excess_n, 1.0, positions + excess_0, &data_temp [ex_excess_0 + ntop + excess_0], &positions_n [0], &data_temp [lda - nbot - ex_excess_n], lda, lda);
					utils::scale (n, alpha_n, &data_temp [lda - 1 - nbot - ex_excess_n - excess_n], lda);
					utils::copy (n, &data_temp [lda - 1], &data_temp [lda - 1 - nbot - ex_excess_n - excess_n], lda, lda);
				} else {
					utils::copy (n, &values_n [0], &data_temp [m - 1 + ntop + ex_excess_0], 1, lda);
				}
					
				if (element_flags & x_solve) {
					TRACE ("Solving in n direction...");
					for (int j = 0; j < m; ++j) {
						utils::diagonal_multiply (n, 1.0, &horizontal_minus_matrix [0], data_in + j, 1.0, &data_temp [ntop + ex_excess_0 + j], 1, m, lda);
						utils::diagonal_solve (n, &horizontal_plus_matrix [0], &data_temp [ntop + ex_excess_0 + j], 1, lda);
					}
					utils::matrix_copy (m, n, &data_temp [ntop + ex_excess_0], data_out, lda);
					
					element_flags &= ~x_solve;
					element_flags |= z_solve;
					
				} else if (element_flags & z_solve) {
					TRACE ("Solving in m direction...");
					utils::matrix_add_scaled (m - 2 + ntop + nbot - excess_0 - excess_n, n, 1.0, data_in + 1 - ntop + excess_0, &data_temp [ex_excess_0 + 1 + excess_0], m, lda);
					utils::interpolate (ex_excess_0, n, m, 1.0, positions, data_in, &positions_0 [0], &data_temp [1], m, lda);
					utils::interpolate (ex_excess_n, n, m, 1.0, positions, data_in, &positions_n [0], &data_temp [lda - 1 - ex_excess_n], m, lda);

					TRACE ("Beginning matrix solve...");
						
					utils::p_block_matrix_solve (messenger_ptr->get_id (), messenger_ptr->get_np (), m - excess_0 - excess_n - ntop - nbot, excess_0 + ex_excess_0 + 2 * ntop, excess_n + ex_excess_n + 2 * nbot, &factorized_matrix [0], &ipiv [0], &data_temp [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, n, lda, sqrt ((int) boundary_matrix.size ()), lda);

					TRACE ("Matrix solve complete.");
					
					for (int i = 0; i < n; ++i) {
						for (int j = 0; j < m; ++j) {
							if (std::isnan (data_temp [i * m + j])) {
								ERROR ("Nan detected.");
								throw 0;
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
					utils::matrix_copy (m, n, &data_temp [ex_excess_0 + ntop], data_out, lda, m);
					
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