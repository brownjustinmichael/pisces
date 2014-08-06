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
#include "source_two_d.hpp"
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
			std::stringstream debug;

			TRACE ("Factorizing...");

			utils::matrix_copy (m, m, &matrix [0], &factorized_matrix [(ex_overlap_0) * (lda + 1)], m, lda);

			/*
				TODO Should we do the matrix copy before the edges?
			*/

			// DEBUG ("ZERO POINT " << &factorized_matrix [0]);

			if (boundary_0) {
				boundary_0->calculate_matrix (timestep, default_matrix + excess_0, &matrix [excess_0], default_matrix, &factorized_matrix [(ex_overlap_0) * (lda + 1) + excess_0], lda);
			}
			if (boundary_n) {
				boundary_n->calculate_matrix (timestep, default_matrix + m - 1 - excess_n, &matrix [m - 1 - excess_n], default_matrix, &factorized_matrix [(ex_overlap_0) * (lda + 1) + m - 1 - excess_n], lda);
			}

			utils::matrix_scale (lda - 2 - excess_0 - ex_overlap_0 - excess_n - ex_overlap_n, lda, timestep, &factorized_matrix [ex_overlap_0 + 1 + excess_0], lda);
			// utils::matrix_scale (lda, lda, timestep, &factorized_matrix [(ex_overlap_0) * (lda + 1)], lda);

			utils::matrix_add_scaled (m - 2 - excess_0 - excess_n, m, 1.0, default_matrix + excess_0 + 1, &factorized_matrix [(ex_overlap_0) * (lda + 1) + excess_0 + 1], m, lda);
			// utils::matrix_add_scaled (m, m, 1.0, default_matrix, &factorized_matrix [(ex_overlap_0) * (lda + 1)], m, lda);

			// for (int j = 0; j < lda; ++j) {
			// 	for (int i = 0; i < lda; ++i) {
			// 		debug << factorized_matrix [i * lda + j] << " ";
			// 	}
			// 	DEBUG (debug.str ());
			// 	debug.str ("");
			// }

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

			// DEBUG ("Solving in m direction..." << &factorized_matrix [0] << " " << &data_temp [0] << " " << &boundary_matrix [0]);

			// for (int j = 0; j < lda; ++j) {
			// 	for (int i = 0; i < ldn; ++i) {
			// 		debug << data_temp [i * lda + j] << " ";
			// 	}
			// 	DEBUG ("RHS " << debug.str ());
			// 	debug.str ("");
			// }

			utils::p_block_matrix_solve (messenger_ptr->get_id (), messenger_ptr->get_np (), inner_m, overlap_0, overlap_n, &factorized_matrix [0], &ipiv [0], &data_temp [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, ldn, lda, sqrt ((int) boundary_matrix.size ()), lda);

			// for (int j = 0; j < lda; ++j) {
			// 	for (int i = 0; i < ldn; ++i) {
			// 		debug << data_temp [i * lda + j] << " ";
			// 	}
			// 	DEBUG ("DONE " << debug.str ());
			// 	debug.str ("");
			// }

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
		bases::solver <datatype> (i_element_flags, i_component_flags), n (i_grid_n.get_n ()), ldn (i_grid_n.get_ld ()), m (i_grid_m.get_n ()), data (i_data), timestep (i_timestep), excess_0 (i_grid_m.get_excess_0 ()), excess_n (i_grid_m.get_excess_n ()), pos_m (&i_grid_m [0]) {
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
		bases::solver <datatype> (i_solver.element_flags, i_solver.component_flags), n (i_solver.grid_ptr (0)->get_n ()),  ldn (i_solver.grid_ptr (0)->get_ld ()),  m (i_solver.grid_ptr (1)->get_n ()), data (i_solver.data_ptr ()), timestep (i_timestep), excess_0 (i_solver.grid_ptr (1)->get_excess_0 ()),  excess_n (i_solver.grid_ptr (1)->get_excess_n ()), pos_m (&((*(i_solver.grid_ptr (1))) [0])) {
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
			std::stringstream debug;
			TRACE ("Executing solve...");
			
			/*
				TODO Add timestep check here?
			*/
			
			// DEBUG ("Solving in n direction");
			
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
			
			// DEBUG ("CHOICE RHS " << rhs_ptr [12 * m + 24] << " " << data [12 * m + 24]);
			//
			// for (int j = 0; j < lda; ++j) {
			// 	for (int i = 0; i < ldn; ++i) {
			// 		debug << data_temp [i * lda + j] << " ";
			// 	}
			// 	DEBUG ("RHS " << debug.str ());
			// 	debug.str ("");
			// }
			//
			// for (int j = 0; j < m; ++j) {
			// 	for (int i = 0; i < ldn; ++i) {
			// 		debug << data [i * m + j] << " ";
			// 	}
			// 	DEBUG ("DATA " << debug.str ());
			// 	debug.str ("");
			// }
			
			for (int j = excess_0; j < m - excess_n; ++j) {
				utils::diagonal_solve (ldn, &factorized_horizontal_matrix [0], &data_temp [ex_overlap_0 + j], 1, lda);
			}
			utils::matrix_copy (m, ldn, &data_temp [ex_overlap_0], data, lda);
			
			for (int i = 0; i < ldn; ++i) {
				for (int j = excess_0 - 1; j >= 0; --j) {
					data [i * m + j] = (data [i * m + j + 2] - data [i * m + j + 1]) / (pos_m [j + 2] - pos_m [j + 1]) * (pos_m [j] - pos_m [j + 1]) + data [i * m + j + 1];
				}
				for (int j = m - excess_n; j < m; ++j) {
					data [i * m + j] = (data [i * m + j - 2] - data [i * m + j - 1]) / (pos_m [j - 2] - pos_m [j - 1]) * (pos_m [j] - pos_m [j - 1]) + data [i * m + j - 1];
				}
			}
			
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
			// DEBUG ("SUP " << &sup [0]);
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
				// DEBUG (data+nbegin);
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
			
			// DEBUG (&sub [nbegin]);
			// DEBUG (&diag [nbegin]);
			// DEBUG (&sup [nbegin]);
			// DEBUG (&supsup [nbegin]);
			// DEBUG (&ipiv [nbegin]);
			// DEBUG (&x [0]);
			// DEBUG (&xipiv [0]);
			// DEBUG (" " << id << " " << np << " " << mm << " " << nbegin);
			utils::p_block_tridiag_solve (id, np, mm, &sub [nbegin], &diag [nbegin], &sup [nbegin], &supsup [nbegin], &ipiv [nbegin], data + nbegin, &x [0], &xipiv [0], &info, ldn, m, m);

			// DEBUG ("DONE");
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
		incompressible_corrector <datatype>::incompressible_corrector (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, utils::messenger* i_messenger_ptr, datatype *i_rhs, datatype* i_data, datatype *i_data_x, datatype *i_data_z, int *i_element_flags, int *i_component_flags, int * i_component_flags_x, int *i_component_flags_z) : 
		bases::solver <datatype> (i_element_flags, i_component_flags),
		n (i_grid_n.get_n ()), 
		ldn (i_grid_n.get_ld ()), 
		m (i_grid_m.get_n ()),
		data (i_data),
		data_x (i_data_x),
		data_z (i_data_z),
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
			component_flags_x = i_component_flags_x;
			component_flags_z = i_component_flags_z;
			
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
				messenger_ptr->send (1, &pos_m [excess_0], id - 1, 0);
			}
			if (id != np - 1) {
				messenger_ptr->recv (1, &ex_pos_m, id + 1, 0);
				messenger_ptr->send (1, &pos_m [m - 2 - excess_n], id + 1, 1);
			}
			if (id != 0) {
				messenger_ptr->recv (1, &ex_pos_0, id - 1, 1);
			}
			data_temp.resize (m * ldn);
			flags = 0x00;
			transform = std::shared_ptr <bases::plan <datatype> > (new fourier::vertical_transform <datatype> (n, m, &data_temp [0], NULL, inverse, i_element_flags, &flags));
			transform_h = std::shared_ptr <bases::plan <datatype> > (new fourier::horizontal_transform <datatype> (n, m, &data_temp [0], NULL, inverse, i_element_flags, &flags));
			z_deriv = std::shared_ptr <bases::plan <datatype >> (new fourier::z_derivative_source <datatype> (grid_n, grid_m, data, data, 1.0, data_z, i_element_flags, i_component_flags_z));
			x_deriv = std::shared_ptr <bases::plan <datatype >> (new fourier::x_derivative_source <datatype> (grid_n, grid_m, data, data, 1.0, data_x, i_element_flags, i_component_flags_z));
		}
		
		template <class datatype>
		incompressible_corrector <datatype>::incompressible_corrector (bases::master_solver <datatype> &i_solver, bases::master_solver <datatype> &i_solver_x, bases::master_solver <datatype> &i_solver_z, utils::messenger* i_messenger_ptr) : incompressible_corrector (*(i_solver.grid_ptr (0)), *(i_solver.grid_ptr (1)), i_messenger_ptr, i_solver.rhs_ptr (spectral_rhs), i_solver.data_ptr (), i_solver_x.data_ptr (), i_solver_z.data_ptr (), i_solver.element_flags, i_solver.component_flags, i_solver_x.component_flags, i_solver_z.component_flags) {}
		
		template <class datatype>
		void incompressible_corrector <datatype>::factorize () {
			TRACE ("Factorizing laplace solver...");

			double scalar = 4.0 * std::acos (-1.0) * std::acos (-1.0) / (pos_n [n - 1] - pos_n [0]) / (pos_n [n - 1] - pos_n [0]);
			int mm = m;
			int nbegin = excess_0;
			if (id != 0) {
				mm -= excess_0 + 1;
			}
			if (id != np - 1) {
				mm -= excess_n + 2;
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
				for (int j = nbegin + 1; j < m - 2 - excess_n; ++j) {
					sub_ptr [i * m + j] = 2.0 / (pos_m [j + 1] - pos_m [j]) / (pos_m [j + 1] - pos_m [j - 1]);
					sup_ptr [i * m + j] = 2.0 / (pos_m [j] - pos_m [j - 1]) / (pos_m [j + 1] - pos_m [j - 1]);
					diag_ptr [i * m + j] = -scalar * (i / 2) * (i / 2) - 2.0 / (pos_m [j + 1] - pos_m [j - 1]) * (1.0 / (pos_m [j + 1] - pos_m [j]) + 1.0 / (pos_m [j] - pos_m [j - 1]));
				}
				if (id != np - 1) {
					sub_ptr [(i + 1) * m - 2 - excess_n] = 2.0 / (ex_pos_m - pos_m [m - 2 - excess_n]) / (ex_pos_m - pos_m [m - 3 - excess_n]);
					sup_ptr [(i + 1) * m - 2 - excess_n] = 2.0 / (pos_m [m - 2 - excess_n] - pos_m [m - 3 - excess_n]) / (ex_pos_m - pos_m [m - 3 - excess_n]);
					diag_ptr [(i + 1) * m - 2 - excess_n] = -scalar * (i / 2) * (i / 2) - 2.0 / (ex_pos_m - pos_m [m - 3 - excess_n]) * (1.0 / (ex_pos_m - pos_m [m - 2 - excess_n]) + 1.0 / (pos_m [m - 2 - excess_n] - pos_m [m - 3 - excess_n]));
				} else {
					sub_ptr [i * m + m - 2] = 2.0 / (pos_m [m - 2 + 1] - pos_m [m - 2]) / (pos_m [m - 2 + 1] - pos_m [m - 2 - 1]);
					sup_ptr [i * m + m - 2] = 2.0 / (pos_m [m - 2] - pos_m [m - 2 - 1]) / (pos_m [m - 2 + 1] - pos_m [m - 2 - 1]);
					diag_ptr [i * m + m - 2] = -scalar * (i / 2) * (i / 2) - 2.0 / (pos_m [m - 2 + 1] - pos_m [m - 2 - 1]) * (1.0 / (pos_m [m - 2 + 1] - pos_m [m - 2]) + 1.0 / (pos_m [m - 2] - pos_m [m - 2 - 1]));
					sup_ptr [(i + 1) * m - 1] = 0.0;
					diag_ptr [(i + 1) * m - 1] = 1.0 / (pos_m [m - 1 - excess_n] - pos_m [m - 2 - excess_n]);
					sub_ptr [(i + 1) * m - 1] = -1.0 / (pos_m [m - 1 - excess_n] - pos_m [m - 2 - excess_n]);
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
			
			std::stringstream debug;
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < ldn; ++i) {
					debug << "SUB " << sub_ptr [i * m + j] << " DIAG " << diag_ptr [i * m + j] << " SUP " << sup_ptr [i * m + j] << " POS " << pos_m [j] << " ";
					
				}
				DEBUG (debug.str ());
				debug.str ("");
			}

			int info;
			
			utils::p_block_tridiag_factorize (id, np, mm, &sub [nbegin], &diag [nbegin], &sup [nbegin], &supsup [nbegin], &ipiv [nbegin], &x [0], &xipiv [0], &info, ldn, m);
		}
		
		template <class datatype>
		void incompressible_corrector <datatype>::execute () {
			std::stringstream debug;
			TRACE ("Solving...");
			int mm = m;
			int nbegin = excess_0;
			if (id != 0) {
				mm -= excess_0 + 1;
			}
			if (id != np - 1) {
				mm -= excess_n + 2;
			}
			
			// if (rhs_ptr) {
			// 	DEBUG (data+nbegin);
			// 	utils::matrix_copy (mm, ldn, rhs_ptr, data + nbegin);
			// }
			
			
			utils::matrix_scale (m, ldn, 0.0, data);
			
			if (z_deriv) {
				z_deriv->execute ();
			}
			
			if (x_deriv) {
				x_deriv->execute ();
			}
			
			
			if (id == 0) {
				utils::scale (ldn, 0.0, data + nbegin, m);
			}
			if (id == np - 1) {
				utils::scale (ldn, 0.0, data + m - 1 - excess_n, m);
			}
			utils::scale (2 * m, 0.0, data);
			
			int info;
			
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < ldn; ++i) {
					debug << data [i * m + j] << " ";
				}
				DEBUG ("RHS " << debug.str ());
				debug.str ("");
			}
			
			utils::p_block_tridiag_solve (id, np, mm, &sub [nbegin], &diag [nbegin], &sup [nbegin], &supsup [nbegin], &ipiv [nbegin], data + nbegin, &x [0], &xipiv [0], &info, ldn, m, m);
			
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < ldn; ++i) {
					debug << data [i * m + j] << " ";
				}
				DEBUG ("DATA " << debug.str ());
				debug.str ("");
			}

			for (int i = 0; i < ldn; ++i) {
				if (id != 0) {
					for (int j = nbegin - excess_0; j >= 0; --j) {
						DEBUG ("EXCESS REGION 0 " << j);
						data [i * m + j] = (data [i * m + j + 2] - data [i * m + j + 1]) / (pos_m [j + 2] - pos_m [j + 1]) * (pos_m [j] - pos_m [j + 1]) + data [i * m + j + 1];
					}
				}
				if (id != np - 1) {
					for (int j = m - excess_n - 1; j < m; ++j) {
						DEBUG ("EXCESS REGION n " << j);
						data [i * m + j] = (data [i * m + j - 2] - data [i * m + j - 1]) / (pos_m [j - 2] - pos_m [j - 1]) * (pos_m [j] - pos_m [j - 1]) + data [i * m + j - 1];
					}
				}
			}
			
			// for (int i = ldn / 3 * 2; i < ldn; ++i) {
			// 	for (int j = 0; j < m; ++j) {
			// 		data [i * m + j] = 0.0;
			// 	}
			// }
			
			utils::matrix_scale (m, ldn, 0.0, &data_temp [0]);
			
			for (int i = 2; i < ldn; i += 2) {
				utils::add_scaled (m, 2.0 * acos (-1.0) / (grid_n [n - 1] - grid_n [0]) * (i / 2), data + i * m, &data_temp [0] + (i + 1) * m);
				utils::add_scaled (m, -2.0 * acos (-1.0) / (grid_n [n - 1] - grid_n [0]) * (i / 2), data + (i + 1) * m, &data_temp [0] + i * m);
			}
			

			
			if (*component_flags_x & transformed_vertical) {
				transform->execute ();
			}
			
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < ldn; ++i) {
					debug << data [i * m + j] << " ";
				}
				DEBUG ("DATA " << debug.str ());
				debug.str ("");
			}
			
			// for (int j = 0; j < m; ++j) {
			// 	for (int i = 0; i < ldn; ++i) {
			// 		debug << data_temp [i * m + j] << " ";
			// 	}
			// 	DEBUG ("GRAD X " << debug.str ());
			// 	debug.str ("");
			// }
			
			// for (int j = 0; j < m; ++j) {
			// 	for (int i = 0; i < ldn; ++i) {
			// 		debug << data_x [i * m + j] << " ";
			// 	}
			// 	DEBUG ("DATA X " << debug.str ());
			// 	debug.str ("");
			// }
			
			utils::matrix_add_scaled (m, ldn, -1.0, &data_temp [0], data_x);
			
			const datatype *pos_m = &grid_m [0];
			datatype *new_data = &data_temp [0];
			
			for (int i = 0; i < ldn; ++i) {
				new_data [i * m] = (data [i * m + 1] - data [i * m]) / (pos_m [1] - pos_m [0]);
				for (int j = 1; j < m - 1; ++j) {
					new_data [i * m + j] = (data [i * m + j + 1] - data [i * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]);
				}
				new_data [(i + 1) * m - 1] = (data [(i + 1) * m - 1] - data [(i + 1) * m - 2]) / (pos_m [m - 1] - pos_m [m - 2]);
			}
			
			if (*component_flags_x & transformed_vertical) {
				transform->execute ();
			}
			
			utils::matrix_add_scaled (m, ldn, -1.0, &data_temp [0], data_z);
			
			utils::scale (2 * m, 0.0, data_z);
			
			utils::matrix_scale (mm, ldn, 0.0, data);
			
			
			if (z_deriv) {
				z_deriv->execute ();
			}
			
			if (x_deriv) {
				x_deriv->execute ();
			}
			
			// for (int j = 0; j < m; ++j) {
			// 	for (int i = 0; i < ldn; ++i) {
			// 		debug << data [i * m + j] << " ";
			// 	}
			// 	DEBUG ("Div U " << debug.str ());
			// 	debug.str ("");
			// }
			
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
		
		template class incompressible_corrector <double>;

	} /* fourier */
} /* two_d */
