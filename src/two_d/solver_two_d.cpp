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
#include "transform_two_d.hpp"
#include "../utils/exceptions.hpp"
#include <sstream>

/*
	TODO boundary class:
	* Hold edge values
	* Communicating boundary can handle overlapping positions
	* - Handle interpolating overlapping zones
	* - Update ntop, nbot
	* Edge boundaries can specify physics of boundary
*/

/*
	TODO Split z, x solvers into separate methods
	* Update base element class
*/

namespace two_d
{
	namespace fourier
	{
		template <class datatype>
		collocation_solver <datatype>::collocation_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, utils::messenger* i_messenger_ptr, datatype& i_timestep, datatype& i_alpha_0, datatype& i_alpha_n, datatype* i_implicit_rhs, datatype *i_explicit_rhs, datatype *i_real_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags) : 
		bases::solver <datatype> (i_element_flags, i_component_flags), n (i_grid_n.get_n ()),  ldn (i_grid_n.get_ld ()),  m (i_grid_m.get_n ()), data (i_data), grid_n (i_grid_n), grid_m (i_grid_m), messenger_ptr (i_messenger_ptr),  timestep (i_timestep),  alpha_0 (i_alpha_0),  alpha_n (i_alpha_n),  positions (&(grid_m [0])), excess_0 (grid_m.get_excess_0 ()),  excess_n (grid_m.get_excess_n ()), default_matrix (grid_m.get_data (0)) {
			TRACE ("Building solver...");
			horizontal_matrix.resize (ldn);
			factorized_horizontal_matrix.resize (ldn);
			matrix.resize (m * m, 0.0);
			values_0.resize (ldn);
			values_n.resize (ldn);
			implicit_rhs_vec = i_implicit_rhs;
			explicit_rhs_vec = i_explicit_rhs;
			real_rhs_vec = i_real_rhs;
			
			ex_excess_0 = 0;
			ntop = 0;
			communicating_boundary <datatype> (messenger_ptr, n, excess_0, ex_excess_0, &positions [0], positions_0, false, ntop);
			// if (messenger_ptr->get_id () - 1 >= 0) {
			// 	messenger_ptr->template send <int> (1, &excess_0, messenger_ptr->get_id () - 1, 0);
			// 	messenger_ptr->template recv <int> (1, &ex_excess_0, messenger_ptr->get_id () - 1, 0);
			// 	positions_0.resize (ex_excess_0);
			// 	messenger_ptr->template send <datatype> (excess_0, positions, messenger_ptr->get_id () - 1, 0);
			// 	messenger_ptr->template recv <datatype> (ex_excess_0, &positions_0 [0], messenger_ptr->get_id () - 1, 0);
			// 	ntop = 1;
			// } else {
			// 	ex_excess_0 = 0;
			// 	ntop = 0;
			// }
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
			ipiv.resize (m); // TODO Should be n - ntop - nbot - excess_0 - excess_n
			data_temp.resize ((m + ex_excess_0 + ex_excess_n + nbot + ntop) * (ldn));
			
			transform = std::shared_ptr <bases::plan <datatype> > (new horizontal_transform <datatype> (n, m, real_rhs_vec, NULL, 0x00, element_flags, &flags));
			/*
				TODO Move this plan to master_solver?
			*/
			TRACE ("Solver built.");
		}
		
		template <class datatype>
		collocation_solver <datatype>::collocation_solver (bases::master_solver <datatype> &i_solver, utils::messenger* i_messenger_ptr, datatype& i_timestep, datatype& i_alpha_0, datatype& i_alpha_n) : 
		bases::solver <datatype> (i_solver.element_flags, i_solver.component_flags), n (i_solver.grid_ptr (0)->get_n ()),  ldn (i_solver.grid_ptr (0)->get_ld ()),  m (i_solver.grid_ptr (1)->get_n ()), data (i_solver.data_ptr ()), grid_n (*(i_solver.grid_ptr (0))), grid_m (*(i_solver.grid_ptr (1))), messenger_ptr (i_messenger_ptr),  timestep (i_timestep),  alpha_0 (i_alpha_0),  alpha_n (i_alpha_n),  positions (&(grid_m [0])), excess_0 (grid_m.get_excess_0 ()),  excess_n (grid_m.get_excess_n ()), default_matrix (grid_m.get_data (0)) {
			TRACE ("Building solver...");
			horizontal_matrix.resize (ldn);
			factorized_horizontal_matrix.resize (ldn);
			matrix.resize (m * m);
			values_0.resize (ldn);
			values_n.resize (ldn);
			implicit_rhs_vec = i_solver.rhs_ptr (implicit_rhs);
			explicit_rhs_vec = i_solver.rhs_ptr (explicit_rhs);
			real_rhs_vec = i_solver.rhs_ptr (real_rhs);
			
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
			ipiv.resize (m); // TODO Should be n - ntop - nbot - excess_0 - excess_n
			data_temp.resize ((m + ex_excess_0 + ex_excess_n + nbot + ntop) * (ldn));
			
			transform = std::shared_ptr <bases::plan <datatype> > (new horizontal_transform <datatype> (n, m, real_rhs_vec, NULL, 0x00, element_flags, &flags));
			/*
				TODO Move this plan to master_solver
			*/
			TRACE ("Solver built.");
		}
		
		template <class datatype>
		void collocation_solver <datatype>::factorize () {
			int info, lda = m + ex_excess_0 + ex_excess_n + nbot + ntop;
			std::stringstream debug;
			TRACE ("Factorizing...");
			
			for (int i = 0; i < ldn; ++i) {
				factorized_horizontal_matrix [i] = 1.0 + timestep * horizontal_matrix [i];
			}
			
			utils::matrix_scale (lda, lda, 0.0, &factorized_matrix [0], lda);
			utils::matrix_copy (m, m, default_matrix, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1)], m, lda);
			
			utils::matrix_add_scaled (m - excess_n - excess_0 - 2, m, timestep, &matrix [0] + excess_0 + 1, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + 1 + excess_0], m, lda);
			if (ntop != 0) {
				utils::matrix_add_scaled (ntop, m, alpha_0 * timestep, &matrix [0] + excess_0, &factorized_matrix [(ntop + ex_excess_0) * lda], m, m + ex_excess_0 + ex_excess_n + ntop + nbot);
				utils::interpolate (ex_excess_0, m, m, timestep, 1.0, positions, &matrix [0], &positions_0 [0], &factorized_matrix [(ntop + ex_excess_0) * lda + ntop], m, lda);
				utils::matrix_add_scaled (ntop, m, alpha_0 * timestep, &matrix [0] + excess_0, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + excess_0], m, m + ex_excess_0 + ex_excess_n + ntop + nbot);
			}
			if (nbot != 0) {
				utils::matrix_add_scaled (nbot, m, alpha_n * timestep, &matrix [0] + m - nbot - excess_n, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + m - nbot - excess_n], m, lda);
				utils::interpolate (ex_excess_n, m, m, timestep, 1.0, positions, &matrix [0], &positions_n [0], &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + m], m, lda);
				utils::matrix_add_scaled (nbot, m, alpha_n * timestep, &matrix [0] + m - nbot - excess_n, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + m + ex_excess_n], m, lda);
			}
			
			utils::p_block_matrix_factorize (messenger_ptr->get_id (), messenger_ptr->get_np (), m - excess_0 - excess_n - ntop - nbot, excess_0 + ex_excess_0 + 2 * ntop, excess_n + ex_excess_n + 2 * nbot, &factorized_matrix [0], &ipiv [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, lda, sqrt ((int) boundary_matrix.size ()));
			
			TRACE ("Done.");
		}
		
		template <class datatype>
		void collocation_solver <datatype>::execute () {
			int info, lda = m + ex_excess_0 + ex_excess_n + nbot + ntop;
			std::stringstream debug;
			TRACE ("Executing solve...");
			
			/*
				TODO Add timestep check here?
			*/
			
			transform->execute ();
			
			utils::scale ((ldn) * lda, 0.0, &data_temp [0]);
			
			if (!(flags & first_run)) {
				utils::copy (ldn, &data [0], &values_0 [0], m);
				utils::copy (ldn, &data [m - 1], &values_n [0], m);
				flags |= first_run;
			}
			
			utils::matrix_add_scaled (m - excess_0 - excess_n, ldn, timestep, &implicit_rhs_vec [excess_0], &data_temp [ex_excess_0 + ntop + excess_0], m, lda);
			utils::matrix_add_scaled (m - excess_0 - excess_n, ldn, timestep, &real_rhs_vec [excess_0], &data_temp [ex_excess_0 + ntop + excess_0], m, lda);
			utils::matrix_add_scaled (m - excess_0 - excess_n, ldn, timestep, &explicit_rhs_vec [excess_0], &data_temp [ex_excess_0 + ntop + excess_0], m, lda);
			
			if (ntop != 0) {
				utils::interpolate (ex_excess_0, ldn, m - excess_0 - excess_n, 1.0, 1.0, positions + excess_0, &data_temp [ex_excess_0 + ntop + excess_0], &positions_0 [0], &data_temp [ntop], lda, lda);
				utils::scale (ldn, alpha_0, &data_temp [0] + ntop + ex_excess_0 + excess_0, lda);
				utils::copy (ldn, &data_temp [ntop + ex_excess_0 + excess_0], &data_temp [0], lda, lda);
			} else {
				utils::copy (ldn, &values_0 [0], &data_temp [0], 1, lda);
			}
			if (nbot != 0) {
				utils::interpolate (ex_excess_n, ldn, m - excess_0 - excess_n, 1.0, 1.0, positions + excess_0, &data_temp [ex_excess_0 + ntop + excess_0], &positions_n [0], &data_temp [lda - nbot - ex_excess_n], lda, lda);
				utils::scale (ldn, alpha_n, &data_temp [0] + lda - 2 * nbot - ex_excess_n - excess_n, lda);
				utils::copy (ldn, &data_temp [0] + lda - 2 * nbot - ex_excess_n - excess_n, &data_temp [0] + lda - nbot, lda, lda);
			} else {
				utils::copy (ldn, &values_n [0], &data_temp [m - 1 + ntop + ex_excess_0], 1, lda);
			}
			
			utils::matrix_add_scaled (m - 2 + ntop + nbot - excess_0 - excess_n, ldn, 1.0, data + 1 - ntop + excess_0, &data_temp [ex_excess_0 + 1 + excess_0], m, lda);
			
			utils::interpolate (ex_excess_0, ldn, m, 1.0, 1.0, positions, data, &positions_0 [0], &data_temp [1], m, lda);
			utils::interpolate (ex_excess_n, ldn, m, 1.0, 1.0, positions, data, &positions_n [0], &data_temp [lda - 1 - ex_excess_n], m, lda);
			
			if (*component_flags & x_solve) {
				TRACE ("Solving in n direction...");
				
				std::vector <datatype> buffer_0 ((excess_0 > ex_excess_0 ? excess_0 + ntop : ex_excess_0 + ntop) * ldn);
				std::vector <datatype> buffer_n ((excess_n > ex_excess_n ? excess_n + nbot : ex_excess_n + nbot) * ldn);
				
				for (int i = 0; i < ldn; ++i) {
					for (int j = 0; j < ex_excess_0 + ntop; ++j) {
						buffer_0 [j * ldn + i] = data_temp [i * lda + j];
					}
				}
				
				if (messenger_ptr->get_id () - 1 >= 0) {
					messenger_ptr->send ((ex_excess_0 + ntop) * ldn, &buffer_0 [0], messenger_ptr->get_id () - 1, 0);
				}
				if (messenger_ptr->get_id () + 1 < messenger_ptr->get_np ()) {
					messenger_ptr->recv ((excess_n + nbot) * ldn, &buffer_n [0], messenger_ptr->get_id () + 1, 0);
				}
				
				for (int i = 0; i < ldn; ++i) {
					for (int j = 0; j < ex_excess_n + nbot; ++j) {
						data_temp [(i + 1) * lda - 2 * nbot - ex_excess_n - excess_n + j] += buffer_n [j * ldn + i];
					}
				}
				for (int i = 0; i < ldn; ++i) {
					for (int j = 0; j < excess_n + nbot; ++j) {
						buffer_n [j * ldn + i] = data_temp [(i + 1) * lda - nbot - ex_excess_n + j];
					}
				}
				
				if (messenger_ptr->get_id () + 1 < messenger_ptr->get_np ()) {
					messenger_ptr->send ((ex_excess_n + nbot) * ldn, &buffer_n [0], messenger_ptr->get_id () + 1, 1);
				}
				if (messenger_ptr->get_id () - 1 >= 0) {
					messenger_ptr->recv ((excess_0 + ntop) * ldn, &buffer_0 [0], messenger_ptr->get_id () - 1, 1);
				}
				
				for (int i = 0; i < ldn; ++i) {
					for (int j = 0; j < ex_excess_0 + ntop; ++j) {
						data_temp [i * lda + ntop + ex_excess_0 + j] += buffer_0 [j * ldn + i];
					}
				}
				
				for (int j = 0; j < m; ++j) {
					utils::diagonal_solve (ldn, &factorized_horizontal_matrix [0], &data_temp [ntop + ex_excess_0 + j], 1, lda);
				}
				utils::matrix_copy (m, ldn, &data_temp [ntop + ex_excess_0], data, lda);
				
			} else if (*component_flags & z_solve) {
				TRACE ("Solving in m direction...");
				
				utils::p_block_matrix_solve (messenger_ptr->get_id (), messenger_ptr->get_np (), m - excess_0 - excess_n - ntop - nbot, excess_0 + ex_excess_0 + 2 * ntop, excess_n + ex_excess_n + 2 * nbot, &factorized_matrix [0], &ipiv [0], &data_temp [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, ldn, lda, sqrt ((int) boundary_matrix.size ()), lda);
				
				TRACE ("Matrix solve complete.");
				
				for (int i = 0; i < ldn; ++i) {
					for (int j = 0; j < m; ++j) {
						if (std::isnan (data_temp [ex_excess_0 + ntop + i * m + j])) {
							FATAL ("Found nan.");
							for (int k = 0; k < m; ++k) {
								printf ("%f ", data_temp [ex_excess_0 + ntop + k * m + j]);
							}
							printf ("\n");
							throw exceptions::nan ();
						}
					}
				}
						
				TRACE ("Updating...");
				utils::matrix_copy (m, ldn, &data_temp [ex_excess_0 + ntop], data, lda, m);
				
				*component_flags |= transformed_vertical;
				
				TRACE ("Solve complete.")
			}
			TRACE ("Execution complete.");
		}
		
		template class collocation_solver <double>;
		
		template <class datatype>
		laplace_solver <datatype>::laplace_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, utils::messenger* i_messenger_ptr, datatype* i_implicit_rhs, datatype *i_explicit_rhs, datatype *i_real_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags) : 
		bases::solver <datatype> (i_element_flags, i_component_flags),
		n (i_grid_n.get_n ()), 
		ldn (i_grid_n.get_ld ()), 
		m (i_grid_m.get_n ()),
		data (i_data),
		grid_n (i_grid_n),
		grid_m (i_grid_m),
		pos_n (&i_grid_n [0]),
		pos_m (&i_grid_m [0]),
		excess_0 (grid_m.get_excess_0 ()), 
		excess_n (grid_m.get_excess_n ()),
		messenger_ptr (i_messenger_ptr) {
			TRACE ("Building laplace solver...");
			sup.resize (m * ldn);
			sub.resize (m * ldn);
			diag.resize (m * ldn);
			implicit_rhs_vec = i_implicit_rhs;
			explicit_rhs_vec = i_explicit_rhs;
			real_rhs_vec = i_real_rhs;
			
			sup_ptr = &sup [0];
			sub_ptr = &sub [0];
			diag_ptr = &diag [0];
			
			supsup.resize (ldn * m);
			ipiv.resize (ldn * m);
			
			if (messenger_ptr->get_id () == 0) {
				x.resize ((m + 4 * messenger_ptr->get_np ()) * 2 * ldn);
				xipiv.resize (2 * messenger_ptr->get_np () * ldn);
			} else {
				x.resize ((m + 4) * 2 * ldn);
			}
			
			id = messenger_ptr->get_id ();
			np = messenger_ptr->get_np ();
			
			if (id != 0) {
				messenger_ptr->send (1, &pos_m [excess_0 + 1], id - 1, 0);
			}
			if (id != np - 1) {
				messenger_ptr->recv (1, &ex_pos_m, id + 1, 0);
				messenger_ptr->send (1, &pos_m [m - 1 - excess_n], id + 1, 1);
			}
			if (id != 0) {
				messenger_ptr->recv (1, &ex_pos_0, id - 1, 1);
			}
			
			transform = std::shared_ptr <bases::plan <datatype> > (new horizontal_transform <datatype> (n, m, real_rhs_vec, NULL, 0x00, i_element_flags, &flags));
		}
		
		template <class datatype>
		laplace_solver <datatype>::laplace_solver (bases::master_solver <datatype> &i_solver, utils::messenger* i_messenger_ptr) : 
		bases::solver <datatype> (i_solver.element_flags, i_solver.component_flags), n ((i_solver.grid_ptr (0))->get_n ()),
		ldn ((i_solver.grid_ptr (0))->get_ld ()),
		m ((i_solver.grid_ptr (1))->get_n ()),
		data (i_solver.data_ptr ()),
		grid_n (*(i_solver.grid_ptr (0))),
		grid_m (*(i_solver.grid_ptr (1))),
		pos_n (&grid_n [0]),
		pos_m (&grid_m [0]),
		excess_0 (grid_m.get_excess_0 ()), 
		excess_n (grid_m.get_excess_n ()),
		messenger_ptr (i_messenger_ptr) {
			TRACE ("Building laplace solver...");
			sup.resize (m * ldn);
			sub.resize (m * ldn);
			diag.resize (m * ldn);
			implicit_rhs_vec = i_solver.rhs_ptr (implicit_rhs);
			explicit_rhs_vec = i_solver.rhs_ptr (explicit_rhs);
			real_rhs_vec = i_solver.rhs_ptr (real_rhs);

			sup_ptr = &sup [0];
			sub_ptr = &sub [0];
			diag_ptr = &diag [0];

			supsup.resize (ldn * m);
			ipiv.resize (ldn * m);

			if (messenger_ptr->get_id () == 0) {
				x.resize ((m + 4 * messenger_ptr->get_np ()) * 2 * ldn);
				xipiv.resize (2 * messenger_ptr->get_np () * ldn);
			} else {
				x.resize ((m + 4) * 2 * ldn);
			}

			id = messenger_ptr->get_id ();
			np = messenger_ptr->get_np ();

			if (id != 0) {
				messenger_ptr->send (1, &pos_m [excess_0 + 1], id - 1, 0);
			}
			if (id != np - 1) {
				messenger_ptr->recv (1, &ex_pos_m, id + 1, 0);
				messenger_ptr->send (1, &pos_m [m - 1 - excess_n], id + 1, 1);
			}
			if (id != 0) {
				messenger_ptr->recv (1, &ex_pos_0, id - 1, 1);
			}

			transform = std::shared_ptr <bases::plan <datatype> > (new horizontal_transform <datatype> (n, m, real_rhs_vec, NULL, 0x00, i_solver.element_flags, &flags));
		}
		
		template <class datatype>
		void laplace_solver <datatype>::factorize () {
			TRACE ("Factorizing laplace solver...");

			double scalar = 4.0 * std::acos (-1.0) / (pos_n [n - 1] - pos_n [0]);
			int mm = m;
			int nbegin = excess_0;
			if (id != 0) {
				mm -= excess_0 + 2;
				nbegin += 1;
			}
			if (id != np - 1) {
				mm -= excess_n + 1;
			}
#pragma omp parallel for
			for (int i = 0; i < ldn; ++i) {
				if (id != 0) {
					sub_ptr [i * m + nbegin] = 2.0 / (pos_m [nbegin + 1] - pos_m [nbegin]) / (pos_m [nbegin + 1] - ex_pos_0);
					sup_ptr [i * m + nbegin] = 2.0 / (pos_m [nbegin] - ex_pos_0) / (pos_m [nbegin + 1] - ex_pos_0);
					diag_ptr [i * m + nbegin] = -scalar * (i / 2) * (i / 2) - 2.0 / (pos_m [nbegin + 1] - ex_pos_0) * (1.0 / (pos_m [nbegin + 1] - pos_m [nbegin]) + 1.0 / (pos_m [nbegin] - ex_pos_0));
				} else {
					sup_ptr [i * m + nbegin] = -1.0 / (pos_m [nbegin + 1] - pos_m [nbegin]);
					diag_ptr [i * m + nbegin] = 1.0 / (pos_m [nbegin + 1] - pos_m [nbegin]);
					sub_ptr [i * m + nbegin] = 0.0;
				}
				for (int j = nbegin + 1; j < m - 1 - excess_n; ++j) {
					sub_ptr [i * m + j] = 2.0 / (pos_m [j + 1] - pos_m [j]) / (pos_m [j + 1] - pos_m [j - 1]);
					sup_ptr [i * m + j] = 2.0 / (pos_m [j] - pos_m [j - 1]) / (pos_m [j + 1] - pos_m [j - 1]);
					diag_ptr [i * m + j] = -scalar * (i / 2) * (i / 2) - 2.0 / (pos_m [j + 1] - pos_m [j - 1]) * (1.0 / (pos_m [j + 1] - pos_m [j]) + 1.0 / (pos_m [j] - pos_m [j - 1]));
				}
				if (id != np - 1) {
					sub_ptr [(i + 1) * m - 1 - excess_n] = 2.0 / (ex_pos_m - pos_m [m - 1 - excess_n]) / (ex_pos_m - pos_m [m - 2 - excess_n]);
					sup_ptr [(i + 1) * m - 1 - excess_n] = 2.0 / (pos_m [m - 1 - excess_n] - pos_m [m - 2 - excess_n]) / (ex_pos_m - pos_m [m - 2 - excess_n]);
					diag_ptr [(i + 1) * m - 1 - excess_n] = -scalar * (i / 2) * (i / 2) - 2.0 / (ex_pos_m - pos_m [m - 2 - excess_n]) * (1.0 / (ex_pos_m - pos_m [m - 1 - excess_n]) + 1.0 / (pos_m [m - 1 - excess_n] - pos_m [m - 2 - excess_n]));
				} else {
					sup_ptr [(i + 1) * m - 1 - excess_n] = 0.0;
					diag_ptr [(i + 1) * m - 1 - excess_n] = 1.0 / (pos_m [m - 1 - excess_n] - pos_m [m - 2 - excess_n]);
					sub_ptr [(i + 1) * m - 1 - excess_n] = -1.0 / (pos_m [m - 1 - excess_n] - pos_m [m - 2 - excess_n]);
				}
			}
			if (id == np - 1) {
				sup_ptr [m - 1] = 0.0;
				sub_ptr [m - 1] = 0.0;
				diag_ptr [m - 1] = 1.0;
				sup_ptr [2 * m - 1] = 0.0;
				sub_ptr [2 * m - 1] = 0.0;
				diag_ptr [2 * m - 1] = 1.0;
			}

			int info;
			
			utils::p_block_tridiag_factorize (id, np, mm, &sub [nbegin], &diag [nbegin], &sup [nbegin], &supsup [nbegin], &ipiv [nbegin], &x [0], &xipiv [0], &info, ldn, m);
		}
		
		template <class datatype>
		void laplace_solver <datatype>::execute () {
			TRACE ("Solving...");
			int mm = m;
			int nbegin = excess_0;
			if (id != 0) {
				mm -= excess_0 + 2;
				nbegin += 1;
			}
			if (id != np - 1) {
				mm -= excess_n + 1;
			}

			transform->execute ();

			utils::matrix_copy (mm, ldn, &explicit_rhs_vec [nbegin], data + nbegin);
			
			utils::matrix_add_scaled (mm, ldn, 1.0, &real_rhs_vec [nbegin], data + nbegin);

			if (id == 0) {
				utils::scale (ldn, 0.0, data + nbegin, m);
			}
			if (id == np - 1) {
				utils::scale (ldn, 0.0, data + m - 1 - excess_n, m);
			}
			utils::scale (2 * m, 0.0, data);

			int info;
			
			utils::p_block_tridiag_solve (id, np, mm, &sub [nbegin], &diag [nbegin], &sup [nbegin], &supsup [nbegin], &ipiv [nbegin], data + nbegin, &x [0], &xipiv [0], &info, ldn, m, m);
			
			for (int i = 0; i < ldn; ++i) {
				for (int j = nbegin - 1; j >= 0; --j) {
					data [i * m + j] = (data [i * m + j + 2] - data [i * m + j + 1]) / (pos_m [j + 2] - pos_m [j + 1]) * (pos_m [j] - pos_m [j + 1]) + data [i * m + j + 1];
				}
				for (int j = m - excess_n; j < m; ++j) {
					data [i * m + j] = (data [i * m + j - 2] - data [i * m + j - 1]) / (pos_m [j - 2] - pos_m [j - 1]) * (pos_m [j] - pos_m [j - 1]) + data [i * m + j - 1];
				}
			}

			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < ldn; ++i) {
					if (std::isnan (data [i * m + j])) {
						FATAL ("Nan in laplace solver.");
						throw exceptions::nan ();
					}
				}
			}

			TRACE ("Solved.");
		}
		
		template class laplace_solver <double>;
	
		template <class datatype>
		divergence_solver <datatype>::divergence_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_x, datatype *i_data_z, int *i_element_flags, int *i_component_flags) : 
		bases::solver <datatype> (i_element_flags, i_component_flags),
		n (i_grid_n.get_n ()), 
		ldn (i_grid_n.get_ld ()), 
		m (i_grid_m.get_n ()),
		pos_m (&(i_grid_m [0])),
		data_x (i_data_x),
		data_z (i_data_z),
		scalar (2.0 * acos (-1.0) / (i_grid_n [n - 1] - i_grid_n [0])),
		grid_n (i_grid_n),
		grid_m (i_grid_m) {}
		
		template <class datatype>
		divergence_solver <datatype>::divergence_solver (bases::master_solver <datatype> &i_solver, datatype *i_data_z) : 
		bases::solver <datatype> (i_solver.element_flags, i_solver.component_flags),
		n (i_solver.grid_ptr (0)->get_n ()), 
		ldn (i_solver.grid_ptr (0)->get_ld ()), 
		m (i_solver.grid_ptr (1)->get_n ()),
		pos_m (&((*(i_solver.grid_ptr (1))) [0])),
		data_x (i_solver.data_ptr ()),
		data_z (i_data_z),
		scalar (2.0 * acos (-1.0) / ((*(i_solver.grid_ptr (0))) [n - 1] - (*(i_solver.grid_ptr (0))) [0])),
		grid_n (*(i_solver.grid_ptr (0))),
		grid_m (*(i_solver.grid_ptr (1))) {}
		
		template <class datatype>
		void divergence_solver <datatype>::factorize () {}
		
		template <class datatype>
		void divergence_solver <datatype>::execute () {
			TRACE ("Solving...");
			utils::scale (m * ldn, 0.0, data_x);
			
			for (int j = 0; j < m; ++j) {
				data_z [j] = 0.0;
				data_z [m + j] = 0.0;
			}
			
			#pragma omp parallel for
			for (int i = 2; i < ldn; i += 2) {
				data_x [i * m] = -(data_z [(i + 1) * m + 1] - data_z [(i + 1) * m]) / (pos_m [1] - pos_m [0]) / scalar / (i / 2);
				data_x [(i + 1) * m] = (data_z [i * m + 1] - data_z [i * m]) / (pos_m [1] - pos_m [0]) / scalar / (i / 2);
				for (int j = 1; j < m - 1; ++j) {
					data_x [i * m + j] = -(data_z [(i + 1) * m + j + 1] - data_z [(i + 1) * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]) / scalar / (i / 2);
					data_x [(i + 1) * m + j] = (data_z [i * m + j + 1] - data_z [i * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]) / scalar / (i / 2);
				}
				data_x [(i + 1) * m - 1] = -(data_z [(i + 2) * m - 1] - data_z [(i + 2) * m - 2]) / (pos_m [m - 1] - pos_m [m - 2]) / scalar / (i / 2);
				data_x [(i + 2) * m - 1] = (data_z [(i + 1) * m - 1] - data_z [(i + 1) * m - 2]) / (pos_m [m - 1] - pos_m [m - 2]) / scalar / (i / 2);
			}
		}
		
		template class divergence_solver <double>;
	} /* fourier */
} /* two_d */
