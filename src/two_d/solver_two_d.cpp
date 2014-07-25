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
		collocation_solver <datatype>::collocation_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, utils::messenger* i_messenger_ptr, datatype& i_timestep, std::shared_ptr <bases::boundary <datatype>> i_boundary_0, std::shared_ptr <bases::boundary <datatype>> i_boundary_n, datatype *i_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags) : 
		bases::solver <datatype> (i_element_flags, i_component_flags), n (i_grid_n.get_n ()), ldn (i_grid_n.get_ld ()), m (i_grid_m.get_n ()), data (i_data), messenger_ptr (i_messenger_ptr), timestep (i_timestep), positions (&(i_grid_m [0])), excess_0 (i_grid_m.get_excess_0 ()), excess_n (i_grid_m.get_excess_n ()), default_matrix (i_grid_m.get_data (0)) {
			TRACE ("Building solver...");
			matrix.resize (m * m, 0.0);
			rhs_ptr = i_rhs;
			
			boundary_0 = i_boundary_0;
			boundary_n = i_boundary_n;
			
			ex_overlap_0 = boundary_0 ? boundary_0->get_ex_overlap () : 0;
			overlap_0 = boundary_0 ? boundary_0->get_overlap () : 0;
			ex_overlap_n = boundary_n ? boundary_n->get_ex_overlap () : 0;
			overlap_n = boundary_n ? boundary_n->get_overlap () : 0;
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
			TRACE ("Solver built.");
		}
		
		template <class datatype>
		collocation_solver <datatype>::collocation_solver (bases::master_solver <datatype> &i_solver, utils::messenger* i_messenger_ptr, datatype& i_timestep, std::shared_ptr <bases::boundary <datatype>> i_boundary_0, std::shared_ptr <bases::boundary <datatype>> i_boundary_n) : 
		bases::solver <datatype> (i_solver.element_flags, i_solver.component_flags), n (i_solver.grid_ptr (0)->get_n ()),  ldn (i_solver.grid_ptr (0)->get_ld ()),  m (i_solver.grid_ptr (1)->get_n ()), data (i_solver.data_ptr ()), messenger_ptr (i_messenger_ptr),  timestep (i_timestep),  positions (&((*(i_solver.grid_ptr (1))) [0])), excess_0 (i_solver.grid_ptr (1)->get_excess_0 ()),  excess_n (i_solver.grid_ptr (1)->get_excess_n ()), default_matrix (i_solver.grid_ptr (1)->get_data (0)) {
			TRACE ("Building solver...");
			matrix.resize (m * m);
			rhs_ptr = i_solver.rhs_ptr (spectral_rhs);
			
			boundary_0 = i_boundary_0;
			boundary_n = i_boundary_n;
			
			ex_overlap_0 = boundary_0 ? boundary_0->get_ex_overlap () : 0;
			overlap_0 = boundary_0 ? boundary_0->get_overlap () : 0;
			ex_overlap_n = boundary_n ? boundary_n->get_ex_overlap () : 0;
			overlap_n = boundary_n ? boundary_n->get_overlap () : 0;
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
			
			/*
				TODO Move this plan to master_solver
			*/
			TRACE ("Solver built.");
		}
		
		template <class datatype>
		void collocation_solver <datatype>::factorize () {
			int info;

			TRACE ("Factorizing...");

			utils::matrix_copy (m, m, &matrix [0], &factorized_matrix [(ex_overlap_0) * (lda + 1)], m, lda);

			/*
				TODO Should we do the matrix copy before the edges?
			*/

			DEBUG ("ZERO POINT " << &factorized_matrix [0]);

			if (boundary_0) {
				boundary_0->calculate_matrix (timestep, default_matrix + excess_0, &matrix [excess_0], &matrix [0], &factorized_matrix [(ex_overlap_0) * (lda + 1) + excess_0], lda);
			}
			if (boundary_n) {
				boundary_n->calculate_matrix (timestep, default_matrix + m - 1 - excess_n, &matrix [m - 1 - excess_n], &matrix [0], &factorized_matrix [(ex_overlap_0) * (lda + 1) + m - 1 - excess_n], lda);
			}

			utils::matrix_scale (lda - 2 - excess_0 - ex_overlap_0 - excess_n - ex_overlap_n, lda, timestep, &factorized_matrix [ex_overlap_0 + 1 + excess_0], lda);
			// utils::matrix_scale (lda, lda, timestep, &factorized_matrix [(ex_overlap_0) * (lda + 1)], lda);

			utils::matrix_add_scaled (m - 2 - excess_0 - excess_n, m, 1.0, default_matrix + excess_0 + 1, &factorized_matrix [(ex_overlap_0) * (lda + 1) + excess_0 + 1], m, lda);
			// utils::matrix_add_scaled (m, m, 1.0, default_matrix, &factorized_matrix [(ex_overlap_0) * (lda + 1)], m, lda);

			utils::p_block_matrix_factorize (messenger_ptr->get_id (), messenger_ptr->get_np (), inner_m, overlap_0, overlap_n, &factorized_matrix [0], &ipiv [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, lda, sqrt ((int) boundary_matrix.size ()));

			TRACE ("Done.");
		}
		
		template <class datatype>
		void collocation_solver <datatype>::execute () {
			int info;
			std::stringstream debug;
			TRACE ("Executing solve...");

			/*
				TODO Add timestep check here?
			*/

			utils::scale ((ldn) * lda, 0.0, &data_temp [0]);

			if (rhs_ptr) {
				utils::matrix_add_scaled (m, ldn, timestep, rhs_ptr, &data_temp [ex_overlap_0], m, lda);
			}

			if (boundary_0) {
				boundary_0->calculate_rhs (data + excess_0, data, &data_temp [0], &data_temp [ex_overlap_0 + excess_0], lda);
			}
			if (boundary_n) {
				boundary_n->calculate_rhs (data + m - 1 - excess_n, data, &data_temp [0], &data_temp [lda - 1 - excess_n - ex_overlap_n], lda);
			}

			utils::matrix_add_scaled (m - 2 - excess_0 - excess_n, ldn, 1.0, data + 1 + excess_0, &data_temp [ex_overlap_0 + 1 + excess_0], m, lda);

			DEBUG ("Solving in m direction..." << &factorized_matrix [0] << " " << &data_temp [0] << " " << &boundary_matrix [0]);

			for (int j = 0; j < lda; ++j) {
				for (int i = 0; i < ldn; ++i) {
					debug << data_temp [i * lda + j] << " ";
				}
				DEBUG ("RHS " << debug.str ());
				debug.str ("");
			}

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
			TRACE ("Execution complete.");
		}
		
		template class collocation_solver <double>;
		
		template <class datatype>
		fourier_solver <datatype>::fourier_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype& i_timestep, std::shared_ptr <bases::boundary <datatype>> i_boundary_0, std::shared_ptr <bases::boundary <datatype>> i_boundary_n, datatype *i_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags) : 
		bases::solver <datatype> (i_element_flags, i_component_flags), n (i_grid_n.get_n ()), ldn (i_grid_n.get_ld ()), m (i_grid_m.get_n ()), data (i_data), timestep (i_timestep), excess_0 (i_grid_m.get_excess_0 ()), excess_n (i_grid_m.get_excess_n ()) {
			TRACE ("Building solver...");
			horizontal_matrix.resize (ldn);
			factorized_horizontal_matrix.resize (ldn);
			rhs_ptr = i_rhs;
			
			boundary_0 = i_boundary_0;
			boundary_n = i_boundary_n;
			
			ex_overlap_0 = boundary_0 ? boundary_0->get_ex_overlap () : 0;
			overlap_0 = boundary_0 ? boundary_0->get_overlap () : 0;
			ex_overlap_n = boundary_n ? boundary_n->get_ex_overlap () : 0;
			overlap_n = boundary_n ? boundary_n->get_overlap () : 0;
			lda = m + ex_overlap_n + ex_overlap_0;
			inner_m = lda - overlap_0 - overlap_n;
			
			data_temp.resize (lda * ldn);
			TRACE ("Solver built.");
		}
		
		template <class datatype>
		fourier_solver <datatype>::fourier_solver (bases::master_solver <datatype> &i_solver, datatype& i_timestep, std::shared_ptr <bases::boundary <datatype>> i_boundary_0, std::shared_ptr <bases::boundary <datatype>> i_boundary_n) : 
		bases::solver <datatype> (i_solver.element_flags, i_solver.component_flags), n (i_solver.grid_ptr (0)->get_n ()),  ldn (i_solver.grid_ptr (0)->get_ld ()),  m (i_solver.grid_ptr (1)->get_n ()), data (i_solver.data_ptr ()), timestep (i_timestep), excess_0 (i_solver.grid_ptr (1)->get_excess_0 ()),  excess_n (i_solver.grid_ptr (1)->get_excess_n ()) {
			TRACE ("Building solver...");
			horizontal_matrix.resize (ldn);
			factorized_horizontal_matrix.resize (ldn);
			rhs_ptr = i_solver.rhs_ptr (spectral_rhs);
			
			boundary_0 = i_boundary_0;
			boundary_n = i_boundary_n;
			
			ex_overlap_0 = boundary_0 ? boundary_0->get_ex_overlap () : 0;
			overlap_0 = boundary_0 ? boundary_0->get_overlap () : 0;
			ex_overlap_n = boundary_n ? boundary_n->get_ex_overlap () : 0;
			overlap_n = boundary_n ? boundary_n->get_overlap () : 0;
			lda = m + ex_overlap_n + ex_overlap_0;
			inner_m = lda - overlap_0 - overlap_n;
			
			data_temp.resize (lda * ldn);
			
			TRACE ("Solver built.");
		}
		
		template <class datatype>
		void fourier_solver <datatype>::factorize () {
			
			TRACE ("Factorizing...");
			
			for (int i = 0; i < ldn; ++i) {
				factorized_horizontal_matrix [i] = 1.0 + timestep * horizontal_matrix [i];
			}
			
			TRACE ("Done.");
		}
		
		template <class datatype>
		void fourier_solver <datatype>::execute () {
			TRACE ("Executing solve...");
			
			/*
				TODO Add timestep check here?
			*/
			
			DEBUG ("Solving in n direction");
			
			utils::scale ((ldn) * lda, 0.0, &data_temp [0]);
			
			if (rhs_ptr) {
				utils::matrix_add_scaled (m, ldn, timestep, rhs_ptr, &data_temp [ex_overlap_0], m, lda);
			}
			
			if (boundary_0) {
				boundary_0->calculate_rhs (data + excess_0, data, &data_temp [0], &data_temp [ex_overlap_0 + excess_0], lda);
			}
			if (boundary_n) {
				boundary_n->calculate_rhs (data + m - 1 - excess_n, data, &data_temp [0], &data_temp [lda - 1 - excess_n - ex_overlap_n], lda);
			}
			
			utils::matrix_add_scaled (m - 2 - excess_0 - excess_n, ldn, 1.0, data + 1 + excess_0, &data_temp [ex_overlap_0 + 1 + excess_0], m, lda);
			
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
			
			TRACE ("Execution complete.");
		}
		
		template class fourier_solver <double>;
		
		template <class datatype>
		laplace_solver <datatype>::laplace_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, utils::messenger* i_messenger_ptr, datatype *i_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags) : 
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
			rhs_ptr = i_rhs;
			
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
			rhs_ptr = i_solver.rhs_ptr (spectral_rhs);

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
			
			if (rhs_ptr) {
				utils::matrix_copy (mm, ldn, rhs_ptr, data + nbegin);
			}
			
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
		horizontal_divergence_solver <datatype>::horizontal_divergence_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_x, datatype *i_data_z, int *i_element_flags, int *i_component_flags) : 
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
		horizontal_divergence_solver <datatype>::horizontal_divergence_solver (bases::master_solver <datatype> &i_solver, datatype *i_data_z) : 
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
		void horizontal_divergence_solver <datatype>::factorize () {}
		
		template <class datatype>
		void horizontal_divergence_solver <datatype>::execute () {
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
		
		template class horizontal_divergence_solver <double>;
		
		// template <class datatype>
		// vertical_divergence_solver <datatype>::vertical_divergence_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, utils::messenger *i_messenger_ptr, datatype* i_data_x, datatype *i_data_z, int *i_element_flags, int *i_component_flags) :
		// bases::solver <datatype> (i_element_flags, i_component_flags),
		// n (i_grid_n.get_n ()),
		// ldn (i_grid_n.get_ld ()),
		// m (i_grid_m.get_n ()),
		// pos_n (&i_grid_n [0]),
		// pos_m (&i_grid_m [0]),
		// data_z (i_data_z),
		// grid_n (i_grid_n),
		// grid_m (i_grid_m),
		// excess_0 (grid_m.get_excess_0 ()),
		// excess_n (grid_m.get_excess_n ()),
		// messenger_ptr (i_messenger_ptr) {
		// 	TRACE ("Building vertical divergence solver...");
		// 	sup.resize (m);
		// 	sub.resize (m);
		// 	diag.resize (m);
		// 	data_x = i_data_x;
		//
		// 	sup_ptr = &sup [0];
		// 	sub_ptr = &sub [0];
		// 	diag_ptr = &diag [0];
		//
		// 	supsup.resize (ldn * m);
		// 	ipiv.resize (ldn * m);
		//
		// 	if (messenger_ptr->get_id () == 0) {
		// 		x.resize ((m + 4 * messenger_ptr->get_np ()) * 2 * ldn);
		// 		xipiv.resize (2 * messenger_ptr->get_np () * ldn);
		// 	} else {
		// 		x.resize ((m + 4) * 2 * ldn);
		// 	}
		//
		// 	id = messenger_ptr->get_id ();
		// 	np = messenger_ptr->get_np ();
		//
		// 	if (id != 0) {
		// 		messenger_ptr->send (1, &pos_m [excess_0 + 1], id - 1, 0);
		// 	}
		// 	if (id != np - 1) {
		// 		messenger_ptr->recv (1, &ex_pos_m, id + 1, 0);
		// 		messenger_ptr->send (1, &pos_m [m - 1 - excess_n], id + 1, 1);
		// 	}
		// 	if (id != 0) {
		// 		messenger_ptr->recv (1, &ex_pos_0, id - 1, 1);
		// 	}
		//
		// 	flags = 0x0;
		// 	transform = std::shared_ptr <bases::plan <datatype> > (new fourier::vertical_transform <datatype> (n, m, data_x, data_z, 0x00, i_element_flags, &flags));
		// }
		//
		// template <class datatype>
		// vertical_divergence_solver <datatype>::vertical_divergence_solver (bases::master_solver <datatype> &i_solver, utils::messenger *i_messenger_ptr, datatype *i_data_x) :
		// bases::solver <datatype> (i_solver.element_flags, i_solver.component_flags),
		// n (i_solver.grid_ptr (0)->get_n ()),
		// ldn (i_solver.grid_ptr (0)->get_ld ()),
		// m (i_solver.grid_ptr (1)->get_n ()),
		// pos_n (&((*(i_solver.grid_ptr (0))) [0])),
		// pos_m (&((*(i_solver.grid_ptr (1))) [0])),
		// data_x (i_data_x),
		// data_z (i_solver.data_ptr ()),
		// grid_n (*(i_solver.grid_ptr (0))),
		// grid_m (*(i_solver.grid_ptr (1))),
		// excess_0 (grid_m.get_excess_0 ()),
		// excess_n (grid_m.get_excess_n ()) {
		// 	TRACE ("Building vertical divergence solver...");
		// 	sup.resize (m);
		// 	sub.resize (m);
		// 	diag.resize (m);
		// 	data_x = i_data_x;
		// 	messenger_ptr = i_messenger_ptr;
		//
		// 	sup_ptr = &sup [0];
		// 	sub_ptr = &sub [0];
		// 	diag_ptr = &diag [0];
		//
		// 	supsup.resize (ldn * m);
		// 	ipiv.resize (ldn * m);
		//
		//
		// 	if (messenger_ptr->get_id () == 0) {
		// 		x.resize ((m + 4 * messenger_ptr->get_np ()) * 2 * ldn);
		// 		xipiv.resize (2 * messenger_ptr->get_np () * ldn);
		// 	} else {
		// 		x.resize ((m + 4) * 2 * ldn);
		// 	}
		//
		// 	id = messenger_ptr->get_id ();
		// 	np = messenger_ptr->get_np ();
		//
		// 	if (id != 0) {
		// 		messenger_ptr->send (1, &pos_m [excess_0 + 1], id - 1, 0);
		// 	}
		// 	if (id != np - 1) {
		// 		messenger_ptr->recv (1, &ex_pos_m, id + 1, 0);
		// 		messenger_ptr->send (1, &pos_m [m - 1 - excess_n], id + 1, 1);
		// 	}
		// 	if (id != 0) {
		// 		messenger_ptr->recv (1, &ex_pos_0, id - 1, 1);
		// 	}
		// 	TRACE ("Done.");
		//
		// 	transform = std::shared_ptr <bases::plan <datatype> > (new fourier::vertical_transform <datatype> (n, m, data_x, data_z, inverse, i_solver.element_flags, &flags));
		// }
		//
		// template <class datatype>
		// void vertical_divergence_solver <datatype>::factorize () {
		// 	TRACE ("Factorizing divergence solver...");
		//
		// 	int mm = m;
		// 	int nbegin = excess_0;
		// 	if (id != 0) {
		// 		mm -= excess_0 + 2;
		// 		nbegin += 1;
		// 	}
		// 	if (id != np - 1) {
		// 		mm -= excess_n + 1;
		// 	}
		// 	if (id != 0) {
		// 		diag_ptr [nbegin] = 1.0 / (pos_m [nbegin] - ex_pos_0);
		// 		sup_ptr [nbegin] = -1.0 / (pos_m [nbegin] - ex_pos_0);
		// 		sub_ptr [nbegin] = 0.0;
		// 	} else {
		// 		sup_ptr [nbegin] = 0.0;
		// 		diag_ptr [nbegin] = 1.0;
		// 		sub_ptr [nbegin] = 0.0;
		// 	}
		// 	for (int j = nbegin + 1; j < m - 1 - excess_n; ++j) {
		// 		diag_ptr [j] = 1.0 / (pos_m [j] - pos_m [j - 1]);
		// 		sup_ptr [j] = -1.0 / (pos_m [j] - pos_m [j - 1]);
		// 		sub_ptr [j] = 0.0;
		// 	}
		// 	if (id != np - 1) {
		// 		diag_ptr [m - 1 - excess_n] = 1.0 / (pos_m [m - 1 - excess_n] - pos_m [m - 2 - excess_n]);
		// 		sup_ptr [m - 1 - excess_n] = -1.0 / (pos_m [m - 1 - excess_n] - pos_m [m - 2 - excess_n]);
		// 		sub_ptr [m - 1 - excess_n] = 0.0;
		// 	} else {
		// 		sup_ptr [m - 1 - excess_n] = 0.0;
		// 		diag_ptr [m - 1 - excess_n] = 1.0;
		// 		sub_ptr [m - 1 - excess_n] = 0.0;
		// 	}
		//
		// 	int info;
		//
		// 	// utils::p_block_tridiag_factorize (id, np, mm, &sub [nbegin], &diag [nbegin], &sup [nbegin], &supsup [nbegin], &ipiv [nbegin], &x [0], &xipiv [0], &info, 1, m);
		// }
		//
		// template <class datatype>
		// void vertical_divergence_solver <datatype>::execute () {
		// 	TRACE ("Solving...");
		// 	int mm = m;
		// 	int nbegin = excess_0;
		// 	if (id != 0) {
		// 		mm -= excess_0 + 2;
		// 		nbegin += 1;
		// 	}
		// 	if (id != np - 1) {
		// 		mm -= excess_n + 1;
		// 	}
		// 	datatype scalar = 2.0 * acos (-1.0) / (pos_n [n - 1] - pos_n [0]);
		// 	flags = transformed_vertical;
		// 	transform->execute ();
		//
		// 	datatype temp;
		// 	if (data_x) {
		// 		// utils::matrix_copy (mm, ldn, data_x, data_z + nbegin);
		// 		for (int i = 2; i < ldn; i += 2) {
		// 			for (int j = 0; j < m; ++j) {
		// 				temp = data_z [i * m + j];
		// 				data_z [i * m + j] = -data_z [(i + 1) * m + j] * scalar * (i / 2);
		// 				data_z [(i + 1) * m + j] = temp * scalar * (i / 2);
		// 			}
		// 		}
		// 	}
		//
		// 	if (id == 0) {
		// 		utils::scale (ldn, 0.0, data_z + nbegin, m);
		// 	}
		// 	if (id == np - 1) {
		// 		utils::scale (ldn, 0.0, data_z + m - 1 - excess_n, m);
		// 	}
		// 	utils::scale (2 * m, 0.0, data_z);
		//
		// 	utils::p_block_direct_tridiag_solve (id, np, mm, &sub [nbegin], &diag [nbegin], &sup [nbegin], &supsup [nbegin], data_z + nbegin, ldn, m);
		//
		//
		// 	for (int i = 0; i < ldn; ++i) {
		// 		for (int j = nbegin - 1; j >= 0; --j) {
		// 			data_z [i * m + j] = (data_z [i * m + j + 2] - data_z [i * m + j + 1]) / (pos_m [j + 2] - pos_m [j + 1]) * (pos_m [j] - pos_m [j + 1]) + data_z [i * m + j + 1];
		// 		}
		// 		for (int j = m - excess_n; j < m; ++j) {
		// 			data_z [i * m + j] = (data_z [i * m + j - 2] - data_z [i * m + j - 1]) / (pos_m [j - 2] - pos_m [j - 1]) * (pos_m [j] - pos_m [j - 1]) + data_z [i * m + j - 1];
		// 		}
		// 	}
		//
		// 	for (int j = 0; j < m; ++j) {
		// 		for (int i = 0; i < ldn; ++i) {
		// 			if (std::isnan (data_z [i * m + j])) {
		// 				FATAL ("Nan in vertical divergence solver.");
		// 				throw exceptions::nan ();
		// 			}
		// 		}
		// 	}
		//
		// 	TRACE ("Solved.");
		// }
		//
		// template class vertical_divergence_solver <double>;
		
		template <class datatype>
		vertical_divergence_solver <datatype>::vertical_divergence_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, utils::messenger* i_messenger_ptr, std::shared_ptr <bases::boundary <datatype>> i_boundary_0, std::shared_ptr <bases::boundary <datatype>> i_boundary_n, datatype *i_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags) :
		bases::solver <datatype> (i_element_flags, i_component_flags), n (i_grid_n.get_n ()), ldn (i_grid_n.get_ld ()), m (i_grid_m.get_n ()), data (i_data), messenger_ptr (i_messenger_ptr), positions (&(i_grid_m [0])), excess_0 (i_grid_m.get_excess_0 ()), excess_n (i_grid_m.get_excess_n ()), default_matrix (i_grid_m.get_data (0)) {
			TRACE ("Building solver...");
			matrix.resize (m * m, 0.0);
			rhs_ptr = i_rhs;

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

			transform = std::shared_ptr <bases::plan <datatype> > (new fourier::vertical_transform <datatype> (n, m, rhs_ptr, &data_temp [ex_overlap_0], inverse, i_element_flags, &flags));
			TRACE ("Solver built.");
		}

		template <class datatype>
		vertical_divergence_solver <datatype>::vertical_divergence_solver (bases::master_solver <datatype> &i_solver, utils::messenger* i_messenger_ptr, std::shared_ptr <bases::boundary <datatype>> i_boundary_0, std::shared_ptr <bases::boundary <datatype>> i_boundary_n, datatype * i_rhs) :
		bases::solver <datatype> (i_solver.element_flags, i_solver.component_flags), n (i_solver.grid_ptr (0)->get_n ()),  ldn (i_solver.grid_ptr (0)->get_ld ()),  m (i_solver.grid_ptr (1)->get_n ()), data (i_solver.data_ptr ()), messenger_ptr (i_messenger_ptr), positions (&((*(i_solver.grid_ptr (1))) [0])), excess_0 (i_solver.grid_ptr (1)->get_excess_0 ()),  excess_n (i_solver.grid_ptr (1)->get_excess_n ()), default_matrix (i_solver.grid_ptr (1)->get_data (0)), deriv_matrix (i_solver.grid_ptr (1)->get_data (1)) {
			TRACE ("Building solver...");
			matrix.resize (m * m);
			rhs_ptr = i_rhs;

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

			transform = std::shared_ptr <bases::plan <datatype> > (new fourier::vertical_transform <datatype> (n, m, rhs_ptr, &data_temp [ex_overlap_0], inverse, i_solver.element_flags, &flags));

			TRACE ("Solver built.");
		}

		template <class datatype>
		void vertical_divergence_solver <datatype>::factorize () {
			int info;

			TRACE ("Factorizing...");

			utils::matrix_scale (m, lda, 0.0, &factorized_matrix [(ex_overlap_0)], lda);

			/*
				TODO Should we do the matrix copy before the edges?
			*/
			std::stringstream debug;
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < m; ++i) {
					debug << factorized_matrix [i * lda + j] << " ";
				}
				DEBUG ("ZERO " << debug.str ());
				debug.str ("");
			}

			DEBUG ("ZERO POINT " << &factorized_matrix [0] << " " << (ex_overlap_0) * (lda + 1) + excess_0);

			if (boundary_0) {
				boundary_0->calculate_matrix (1.0, default_matrix + excess_0, deriv_matrix + excess_0, deriv_matrix, &factorized_matrix [(ex_overlap_0) * (lda + 1) + excess_0], lda, true);
			}
			if (boundary_n) {
				boundary_n->calculate_matrix (1.0, default_matrix + m - 1 - excess_n, deriv_matrix + m - 1 - excess_n, deriv_matrix, &factorized_matrix [(ex_overlap_0) * (lda + 1) + m - 1 - excess_n], lda, true);
			}

			// utils::matrix_scale (lda - 2 - excess_0 - ex_overlap_0 - excess_n - ex_overlap_n, lda, timestep, &factorized_matrix [ex_overlap_0 + 1 + excess_0], lda);

			utils::matrix_add_scaled (m - 2 - excess_0 - excess_n, m, 1.0, deriv_matrix + excess_0 + 1, &factorized_matrix [(ex_overlap_0) * (lda + 1) + excess_0 + 1], m, lda);
			
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < m; ++i) {
					debug << deriv_matrix [i * m + j] << " ";
				}
				DEBUG ("DERIV " << debug.str ());
				debug.str ("");
			}
			for (int j = 0; j < lda; ++j) {
				for (int i = 0; i < lda; ++i) {
					debug << factorized_matrix [i * lda + j] << " ";
				}
				DEBUG (debug.str ());
				debug.str ("");
			}

			utils::p_block_matrix_factorize (messenger_ptr->get_id (), messenger_ptr->get_np (), inner_m, overlap_0, overlap_n, &factorized_matrix [0], &ipiv [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, lda, sqrt ((int) boundary_matrix.size ()));

			TRACE ("Done.");
		}

		template <class datatype>
		void vertical_divergence_solver <datatype>::execute () {
			int info;
			std::stringstream debug;
			TRACE ("Executing solve...");

			/*
				TODO Add timestep check here?
			*/

			utils::scale ((ldn) * lda, 0.0, &data_temp [0]);

			transform->execute ();

			if (boundary_0) {
				boundary_0->calculate_rhs (NULL, NULL, &data_temp [0], &data_temp [ex_overlap_0 + excess_0], lda);
			}
			if (boundary_n) {
				boundary_n->calculate_rhs (NULL, NULL, &data_temp [0], &data_temp [lda - 1 - excess_n - ex_overlap_n], lda);
			}

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
			TRACE ("Execution complete.");
		}

		template class vertical_divergence_solver <double>;
	} /* fourier */
} /* two_d */
