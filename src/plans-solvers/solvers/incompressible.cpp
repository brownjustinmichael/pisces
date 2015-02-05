/*!**********************************************************************
 * \file solver_two_d.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-14.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "incompressible.hpp"

#include <cmath>
#include <sstream>

#include "linalg/utils.hpp"
#include "linalg/linalg.hpp"
#include "linalg/interpolate.hpp"
#include "linalg/exceptions.hpp"
#include "linalg-block/banded.hpp"

#include "plans-transforms/transform.hpp"

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

namespace plans
{
	template <class datatype>
	incompressible_corrector <datatype>::incompressible_corrector (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, mpi::messenger* i_messenger_ptr, std::shared_ptr <plans::boundary <datatype>> i_boundary_0, std::shared_ptr <plans::boundary <datatype>> i_boundary_n, datatype *i_rhs, datatype* i_data, datatype *i_data_x, datatype *i_data_z, int *i_element_flags, int *i_component_flags, int * i_component_flags_x, int *i_component_flags_z) :
	plans::solver <datatype> (i_element_flags, i_component_flags),
	n (i_grid_n.get_n ()),
	ldn (i_grid_n.get_ld ()),
	m (i_grid_m.get_n ()),
	data (i_data),
	data_x (i_data_x),
	data_z (i_data_z),
	grid_n (i_grid_n),
	grid_m (i_grid_m),
	pos_n (&grid_n [0]),
	excess_0 (grid_m.get_excess_0 ()),
	excess_n (grid_m.get_excess_n ()),
	messenger_ptr (i_messenger_ptr) {
		TRACE ("Building laplace solver...");
		rhs_ptr = i_rhs;
		component_flags_x = i_component_flags_x;
		component_flags_z = i_component_flags_z;
		
		boundary_0 = i_boundary_0;
		boundary_n = i_boundary_n;

		ipiv.resize (ldn * (m + 2));

		if (messenger_ptr->get_id () == 0) {
			x.resize (4 * (3) * (3) * ldn * messenger_ptr->get_np () * messenger_ptr->get_np ());
			xipiv.resize (2 * (3) * messenger_ptr->get_np () * ldn);
		} else {
			x.resize (4 * (3) * (3) * ldn);
			xipiv.resize (1);
		}

		id = messenger_ptr->get_id ();
		np = messenger_ptr->get_np ();
		
		positions.resize (m + 6);
		mid_positions.resize (m + 6);
		diff_pos = &positions [3];
		
		for (int j = 0; j < m - 0; ++j) {
			diff_pos [j] = grid_m [j];
		}
		
		ntop = 0;
		nbot = 0;
		if (id != 0) {
			ntop = 1;
			messenger_ptr->send (3, &diff_pos [excess_0 + 1], id - 1, 0);
		}
		if (id != np - 1) {
			nbot = 2;
			messenger_ptr->recv (3, &diff_pos [m - excess_n], id + 1, 0);
			messenger_ptr->send (3, &diff_pos [m - 4 - excess_n], id + 1, 1);
		} else {
			diff_pos [m] = diff_pos [m - 1] + (diff_pos [m - 1] - diff_pos [m - 2]);
		}
		if (id != 0) {
			messenger_ptr->recv (3, &diff_pos [excess_0 - 3], id - 1, 1);
		} else {
			diff_pos [-1] = diff_pos [0] - (diff_pos [1] - diff_pos [0]);
		}
		
		diff_midpos = &mid_positions [3];
		for (int j = -3; j < m + 3; ++j) {
			diff_midpos [j] = (diff_pos [j] + diff_pos [j + 1]) / 2.0;
		}
		if (id == 0) {
			diff_midpos [-1] = 2.0 * diff_pos [0] - (diff_pos [0] + diff_pos [1]) / 2.0;
		}
		if (id == np - 1) {
			diff_midpos [m - 1] = 2.0 * diff_pos [m - 1] - (diff_pos [m - 1] + diff_pos [m - 2]) / 2.0;
			diff_midpos [m] = diff_pos [m - 1] + (diff_pos [m - 1] - diff_midpos [m - 3]);
		}
		
		for (int j = m + 2; j >= -3; --j) {
			diff_midpos [j] -= diff_midpos [j - 1];
			diff_pos [j] -= diff_pos [j - 1];
		}
		
		diff_midpos += excess_0;
		diff_pos += excess_0;
		
		int kl = 2;
		int ku = 1;
		int lda = 2 * kl + ku + 1;
		
		matrix.resize ((lda) * (m + 2 + kl + ku) * ldn);
		bufferl.resize (np * (m + 2) * kl * ldn);
		bufferr.resize (np * (m + 2) * ku * ldn);

		data_temp.resize ((m + 2) * ldn);
		
		count = 0;
	}

	template <class datatype>
	void incompressible_corrector <datatype>::factorize () {
		std::stringstream debug;
		int info;
		TRACE ("Factorizing laplace solver...");
		if (count == 0) {

			double scalar = 4.0 * std::acos (-1.0) * std::acos (-1.0) / (pos_n [n - 1] - pos_n [0]) / (pos_n [n - 1] - pos_n [0]);
			datatype *matrix_ptr;
			int kl = 2;
			int ku = 1;
			int lda = 2 * kl + ku + 1;
			
			for (int i = 0; i < ldn; ++i) {
				matrix_ptr = &matrix [(i) * (m + 2 + kl + ku) * lda + kl + ku + (kl + ku + excess_0) * lda];
				for (int j = 0; j < m + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0); ++j) {
					// j is the gridpoint location of the solve on the velocity grid (I think)
					matrix_ptr [(j - 2) * lda + 2] = 1.0 / diff_midpos [j - 1] / diff_pos [j + 1];
					matrix_ptr [(j - 1) * lda + 1] = -1.0 / diff_midpos [j - 1] / diff_pos [j + 1] - scalar * (i / 2) * (i / 2) / 2.0;
					matrix_ptr [(j) * lda] = -1.0 / diff_midpos [j + 1] / diff_pos [j + 1] - scalar * (i / 2) * (i / 2) / 2.0;
					matrix_ptr [(j + 1) * lda - 1] = 1.0 / diff_midpos [j + 1] / diff_pos [j + 1];
				}
				if (id == 0) {
					matrix_ptr [-3 * lda + 2] = 0.0;
					matrix_ptr [-2 * lda + 1] = 0.0;
					matrix_ptr [-1 * lda] = 1.0;
					matrix_ptr [-1] = -1.0;

					matrix_ptr [-2 * lda + 2] = 0.0;
					matrix_ptr [-1 * lda + 1] = 1.0 / diff_midpos [0] / (diff_pos [1]);
					matrix_ptr [0] = -1.0 / diff_midpos [1] / diff_pos [1] - 1.0 / diff_midpos [0] / diff_pos [1] - scalar * (i / 2) * (i / 2);
					matrix_ptr [1 * lda - 1] = 1.0 / diff_midpos [1] / diff_pos [1];
				}
				if (id == np - 1) {
					int j = m + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0);
					matrix_ptr [(j - 2) * lda + 2] = -1.0;
					matrix_ptr [(j - 1) * lda + 1] = 1.0;
					matrix_ptr [j * lda] = 1.0e-10;
					matrix_ptr [(j + 1) * lda - 1] = 0.0;
				}
			}

			linalg::p_block_banded_factorize (id, np, m + (nbot == 0 ? 1 : -nbot - excess_n - 1) + (id == 0 ? 1: -excess_0 - ntop), kl, ku, &matrix [(id == 0 ? 0 : 1 + excess_0) * lda], &ipiv [0], &x [0], &xipiv [0], &bufferl [0], &bufferr [0], &info, ldn, lda, m + 2 + kl + ku);
		}

	}

	template <class datatype>
	void incompressible_corrector <datatype>::execute () {
		std::stringstream debug;
		int info;
		TRACE ("Solving...");
		bool retransform = false;
		int kl = 2;
		int ku = 1;
		int lda = 2 * kl + ku + 1;
		
		DEBUG ("HERE")
		
		linalg::scale ((m + 2) * ldn, 0.0, &data_temp [0]);
		
		// std::shared_ptr <plans::plan <datatype> > transform_x = std::shared_ptr <plans::plan <datatype> > (new plans::vertical_transform <datatype> (n, m, data_x, NULL, 0x00, element_flags, component_flags_x));
		// std::shared_ptr <plans::plan <datatype> > transform_z = std::shared_ptr <plans::plan <datatype> > (new plans::vertical_transform <datatype> (n, m, data_z, NULL, 0x00, element_flags, component_flags_z));
		
		if (*component_flags_x & transformed_vertical) {
			// transform_x->execute ();
			// transform_z->execute ();
			// retransform = true;
		}
		if (!(*component_flags_x & transformed_vertical)) {
			datatype scalar = acos (-1.0) * 2.0 / (pos_n [n - 1] - pos_n [0]);
			datatype *data_ptr = &data_temp [1];
			int mmax = m + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0);

			datatype *ndata_z = &data_z [excess_0], *ndata_x = &data_x [excess_0];

			data_ptr += excess_0;
			
			linalg::scale (ldn * (m + 2), 0.0, &data_temp [0]);
			
			for (int i = 2; i < ldn; i += 2) {
				for (int j = -1; j < mmax + 1; ++j) {
					data_ptr [i * (m + 2) + j] = -scalar * (i / 2) * ndata_x [(i + 1) * m + j];
					data_ptr [(i + 1) * (m + 2) + j] = scalar * (i / 2) * ndata_x [i * m + j];
				}
			}
			
			// For every point, we need to add the vertical derivative of the vertical component to the right hand side 
			for (int j = (id == 0 ? 1 : 0); j < mmax; ++j) {
				linalg::add_scaled (ldn, 1.0 / diff_pos [j + 1], ndata_z + j + 1, data_ptr + j, m, m + 2);
				linalg::add_scaled (ldn, -1.0 / diff_pos [j + 1], ndata_z + j - 1, data_ptr + j, m, m + 2);
			}
			
			// The vertical boundaries are special cases for the vertical derivatives
			for (int i = 0; i < ldn; ++i) {
				if (id == 0) {
					data_ptr [i * (m + 2) - 1] = 0.0;
					data_ptr [i * (m + 2)] += (ndata_z [i * m + 1] - ndata_z [i * m]) / (diff_pos [1]);
				}
				if (id == np - 1) {
					data_ptr [i * (m + 2) + (m + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0))] = 0.0;
				}
			}
			
			// Solve the matrix
			linalg::p_block_banded_solve (id, np, m + (nbot == 0 ? 1 : -nbot - excess_n - 1) + (id == 0 ? 1: -excess_0 - ntop), kl, ku, &matrix [(id == 0 ? 0 : 1 + excess_0) * 6], &ipiv [0], &data_temp [(id == 0 ? 0 : 1 + excess_0)], &x [0], &xipiv [0], &bufferl [0], &bufferr [0], &info, ldn, lda, m + 2 + kl + ku, m + 2);
			
			// We have chosen to be in the frame with no net vertical velocity (we can always choose this)
			linalg::scale (2 * (m + 2), 0.0, &data_temp [0]);
			
			// Communicate the resulting boundary conditions
			if (boundary_n) {
				boundary_n->send (&data_temp [m - excess_n - 1], m + 2, 1);
			}
			if (boundary_0) {
				boundary_0->receive (&data_ptr [-1], m + 2, 1, 0.0);
				boundary_0->send (&data_ptr [0], m + 2, 1);
			}
			if (boundary_n) {
				boundary_n->receive (&data_temp [m - excess_n], m + 2, 1, 0.0);
			}
			
			// For plotting purposes, we can calculate the resulting pressure
			linalg::matrix_copy (m, ldn, data_ptr, data, m + 2);
			linalg::matrix_add_scaled (m, ldn, 1.0, data_ptr + 1, data, m + 2);
			linalg::scale (m * ldn, 0.5, data);
			
			// The vertical component incompressible correction
			for (int j = 0; j < mmax; ++j) {
				linalg::add_scaled (ldn - 2, -1.0 / diff_midpos [j], data_ptr + 2 * (m + 2) + j, ndata_z + 2 * m + j, m + 2, m);
				linalg::add_scaled (ldn - 2, 1.0 / diff_midpos [j], data_ptr + 2 * (m + 2) + j - 1, ndata_z + 2 * m + j, m + 2, m);
			}
			
			// The horizontal componenet incompressible correction
			for (int i = 2; i < ldn; i += 2) {
				linalg::add_scaled (mmax, scalar * (i / 2) / 2.0, data_ptr + (i + 1) * (m + 2) - 1, ndata_x + i * m);
				linalg::add_scaled (mmax, scalar * (i / 2) / 2.0, data_ptr + (i + 1) * (m + 2), ndata_x + i * m);
				linalg::add_scaled (mmax, -scalar * (i / 2) / 2.0, data_ptr + i * (m + 2) - 1, ndata_x + (i + 1) * m);
				linalg::add_scaled (mmax, -scalar * (i / 2) / 2.0, data_ptr + i * (m + 2), ndata_x + (i + 1) * m);
			}
			
			// Since only one boundary was solved, get the other and the overlap region from the adjacent element
			if (boundary_0) {
				// This should be 1 + ex_excess_0
				boundary_0->send (&data_z [excess_0], m, 1 + excess_0);
			}
			if (boundary_n) {
				boundary_n->receive (&data_z [m - 1 - excess_n], m, 1 + excess_n, 0.0);
				// This should be m - 1 - excess_n - ex_excess_n and ex_excess_n
				boundary_n->send (&data_z [m - 1 - excess_n - excess_n], m, excess_n);
			}
			if (boundary_0) {
				boundary_0->receive (&data_z [0], m, excess_0, 0.0);
			}

			if (boundary_0) {
				// This should be 1 + ex_excess_0
				boundary_0->send (&data_x [excess_0], m, 1 + excess_0);
			}
			if (boundary_n) {
				boundary_n->receive (&data_x [m - 1 - excess_n], m, 1 + excess_n, 0.0);
				// This should be m - 1 - excess_n - ex_excess_n and ex_excess_n
				boundary_n->send (&data_x [m - 1 - excess_n - excess_n], m, excess_n);
			}
			if (boundary_0) {
				boundary_0->receive (&data_x [0], m, excess_0, 0.0);
			}
			
			linalg::scale (2 * m, 0.0, data_z);
			linalg::scale (2 * m, 0.0, data_x);

#ifdef CHECKNAN
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < ldn; ++i) {
					if (std::isnan (data [i * m + j])) {
						FATAL ("Nan in laplace solver.");
						throw linalg::exceptions::nan ();
					}
				}
			}
#endif
		} else {
			// FATAL ("SHOULDN'T BE HERE" << *component_flags_x);
			// throw 0;
		}
		
		// if (retransform) {
		// 	transform_x->execute ();
		// 	transform_z->execute ();
		// }

		TRACE ("Solved");
	}

	template class incompressible_corrector <double>;
} /* plans */
