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
		collocation_solver <datatype>::collocation_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, utils::messenger* i_messenger_ptr, datatype& i_timestep, std::shared_ptr <bases::boundary <datatype>> i_boundary_0, std::shared_ptr <bases::boundary <datatype>> i_boundary_n, datatype* i_implicit_rhs, datatype *i_explicit_rhs, datatype *i_real_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags) : 
		bases::solver <datatype> (i_element_flags, i_component_flags), n (i_grid_n.get_n ()), ldn (i_grid_n.get_ld ()), m (i_grid_m.get_n ()), data (i_data), grid_n (i_grid_n), grid_m (i_grid_m), messenger_ptr (i_messenger_ptr), timestep (i_timestep), positions (&(grid_m [0])), excess_0 (grid_m.get_excess_0 ()), excess_n (grid_m.get_excess_n ()), default_matrix (grid_m.get_data (0)) {
			TRACE ("Building solver...");
			horizontal_matrix.resize (ldn);
			factorized_horizontal_matrix.resize (ldn);
			matrix.resize (m * m, 0.0);
			implicit_rhs_vec = i_implicit_rhs;
			explicit_rhs_vec = i_explicit_rhs;
			real_rhs_vec = i_real_rhs;
			
			boundary_0_ptr = i_boundary_0;
			boundary_0 = &*boundary_0_ptr;
			boundary_n_ptr = i_boundary_n;
			boundary_n = &*boundary_n_ptr;
			
			ex_overlap_0 = boundary_0->get_ex_overlap ();
			overlap_0 = boundary_0->get_overlap ();
			ex_overlap_n = boundary_n->get_ex_overlap ();
			overlap_n = boundary_n->get_overlap ();
			lda = m + ex_overlap_n + ex_overlap_0;
			inner_m = lda - overlap_0 - overlap_n;
			int ns0 = overlap_0;
			
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
				boundary_matrix.resize ((overlap_0 + overlap_n) * (overlap_0 + overlap_n));
			}
			
			factorized_matrix.resize (lda * lda, 0.0);
			ipiv.resize (m); // TODO Should be n - ntop - nbot - excess_0 - excess_n
			data_temp.resize (lda * ldn);
			
			transform = std::shared_ptr <bases::plan <datatype> > (new horizontal_transform <datatype> (n, m, real_rhs_vec, NULL, 0x00, element_flags, &flags));
			/*
				TODO Move this plan to master_solver?
			*/
			TRACE ("Solver built.");
		}
		
		template <class datatype>
		collocation_solver <datatype>::collocation_solver (bases::master_solver <datatype> &i_solver, utils::messenger* i_messenger_ptr, datatype& i_timestep, std::shared_ptr <bases::boundary <datatype>> i_boundary_0, std::shared_ptr <bases::boundary <datatype>> i_boundary_n) : 
		bases::solver <datatype> (i_solver.element_flags, i_solver.component_flags), n (i_solver.grid_ptr (0)->get_n ()),  ldn (i_solver.grid_ptr (0)->get_ld ()),  m (i_solver.grid_ptr (1)->get_n ()), data (i_solver.data_ptr ()), grid_n (*(i_solver.grid_ptr (0))), grid_m (*(i_solver.grid_ptr (1))), messenger_ptr (i_messenger_ptr),  timestep (i_timestep),  positions (&(grid_m [0])), excess_0 (grid_m.get_excess_0 ()),  excess_n (grid_m.get_excess_n ()), default_matrix (grid_m.get_data (0)) {
			TRACE ("Building solver...");
			horizontal_matrix.resize (ldn);
			factorized_horizontal_matrix.resize (ldn);
			matrix.resize (m * m);
			implicit_rhs_vec = i_solver.rhs_ptr (implicit_rhs);
			explicit_rhs_vec = i_solver.rhs_ptr (explicit_rhs);
			real_rhs_vec = i_solver.rhs_ptr (real_rhs);
			
			boundary_0 = i_boundary_0;
			boundary_n = i_boundary_n;
			
			ex_overlap_0 = boundary_0->get_ex_overlap ();
			overlap_0 = boundary_0->get_overlap ();
			ex_overlap_n = boundary_n->get_ex_overlap ();
			overlap_n = boundary_n->get_overlap ();
			lda = m + ex_overlap_n + ex_overlap_0;
			inner_m = lda - overlap_0 - overlap_n;
			
			int ns0 = overlap_0;
			
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
				boundary_matrix.resize ((overlap_0 + overlap_n) * (overlap_0 + overlap_n));
			}
			
			factorized_matrix.resize (lda * lda, 0.0);
			ipiv.resize (m); // TODO Should be n - ntop - nbot - excess_0 - excess_n
			data_temp.resize (lda * ldn);
			
			transform = std::shared_ptr <bases::plan <datatype> > (new horizontal_transform <datatype> (n, m, real_rhs_vec, NULL, 0x00, element_flags, &flags));
			/*
				TODO Move this plan to master_solver
			*/
			TRACE ("Solver built.");
		}
		
		template <class datatype>
		void collocation_solver <datatype>::factorize () {
			int info;
			std::stringstream debug;
			
			TRACE ("Factorizing...");
			
			for (int i = 0; i < ldn; ++i) {
				factorized_horizontal_matrix [i] = 1.0 + timestep * horizontal_matrix [i];
			}
			
			utils::matrix_copy (m, m, &matrix [0], &factorized_matrix [(ex_overlap_0) * (lda + 1)], m, lda);
			if (boundary_0) {
				boundary_0->calculate_matrix (&matrix [excess_0], &matrix [0], &factorized_matrix [(ex_overlap_0) * lda], lda);
			}
			if (boundary_n) {
				boundary_n->calculate_matrix (&matrix [m - ex_overlap_n], &matrix [0], &factorized_matrix [(ex_overlap_0) * (lda + 1) + m - ex_overlap_n], lda);
			}
			
			utils::matrix_scale (lda, lda, timestep, &factorized_matrix [0], lda);
			
			utils::matrix_add_scaled (m, m, 1.0, default_matrix, &factorized_matrix [(ex_overlap_0) * (lda + 1)], m, lda);
			
			utils::p_block_matrix_factorize (messenger_ptr->get_id (), messenger_ptr->get_np (), inner_m, overlap_0, overlap_n, &factorized_matrix [0], &ipiv [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, lda, sqrt ((int) boundary_matrix.size ()));
			
			TRACE ("Done.");
		}
		
		template <class datatype>
		void collocation_solver <datatype>::execute () {
			int info;
			TRACE ("Executing solve...");
			
			/*
				TODO Add timestep check here?
			*/
			
			transform->execute ();
			
			utils::scale ((ldn) * lda, 0.0, &data_temp [0]);
			
			utils::matrix_add_scaled (m, ldn, timestep, &implicit_rhs_vec [0], &data_temp [ex_overlap_0], m, lda);
			utils::matrix_add_scaled (m, ldn, timestep, &real_rhs_vec [0], &data_temp [ex_overlap_0], m, lda);
			utils::matrix_add_scaled (m, ldn, timestep, &explicit_rhs_vec [0], &data_temp [ex_overlap_0], m, lda);
			
			/*
				TODO Once the transform is moved to master_solver, there's no need for these to be separate
			*/
			
			if (boundary_0) {
				boundary_0->calculate_rhs (data, &data_temp [0], &data_temp [0], lda);
			}
			if (boundary_n) {
				boundary_n->calculate_rhs (data, &data_temp [0], &data_temp [lda - 1 - excess_n - ex_overlap_n], lda);
			}
			
			utils::matrix_add_scaled (m - 2 - excess_0 - excess_n, ldn, 1.0, data + 1 + excess_0, &data_temp [ex_overlap_0 + 1 + excess_0], m, lda);
			
			if (*component_flags & x_solve) {
				TRACE ("Solving in n direction...");
				
				if (boundary_0) {
					boundary_0->send (&data_temp [0], lda);
				}
				if (boundary_n) {
					boundary_n->receive (&data_temp [lda - overlap_n], lda);
					boundary_n->send (&data_temp [lda - ex_overlap_n], lda);
				}
				if (boundary_0) {
					boundary_0->receive (&data_temp [ex_overlap_0], lda);
				}
				
				for (int j = 0; j < m; ++j) {
					utils::diagonal_solve (ldn, &factorized_horizontal_matrix [0], &data_temp [ex_overlap_0 + j], 1, lda);
				}
				utils::matrix_copy (m, ldn, &data_temp [ex_overlap_0], data, lda);
				
			} else if (*component_flags & z_solve) {
				TRACE ("Solving in m direction...");
				
				utils::p_block_matrix_solve (messenger_ptr->get_id (), messenger_ptr->get_np (), inner_m, overlap_0, overlap_n, &factorized_matrix [0], &ipiv [0], &data_temp [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, ldn, lda, sqrt ((int) boundary_matrix.size ()), lda);
				
				TRACE ("Matrix solve complete.");
				
				for (int i = 0; i < ldn; ++i) {
					for (int j = 0; j < m; ++j) {
						if (std::isnan (data_temp [ex_overlap_0 + i * m + j])) {
							FATAL ("Found nan.");
							for (int k = 0; k < m; ++k) {
								printf ("%f ", data_temp [ex_overlap_0 + k * m + j]);
							}
							printf ("\n");
							throw exceptions::nan ();
						}
					}
				}
						
				TRACE ("Updating...");
				utils::matrix_copy (m, ldn, &data_temp [ex_overlap_0], data, lda, m);
				
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
