/*!**********************************************************************
 * \file pseudo.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-14.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "pseudo.hpp"

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
	* - Update ntop, nbot
	* Edge boundaries can specify physics of boundary
*/

namespace plans
{
	namespace solvers
	{
		template <class datatype>
		pseudo_incompressible <datatype>::pseudo_incompressible (mpi::messenger* i_messenger_ptr, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_0, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_n, grids::variable <datatype> &i_data, grids::variable <datatype> &i_data_out, grids::variable <datatype> &i_rhs, grids::variable <datatype> &i_data_x, grids::variable <datatype> &i_data_z, datatype *i_density, datatype *i_pressure, datatype i_gamma) : 
		solver <datatype> (i_data, i_data_out, this->get_state_in (), this->get_state ()), 
		n (i_data.get_grid (0).get_n ()), 
		ldn (i_data.get_grid (0).get_ld ()), 
		m (i_data.get_grid (1).get_n ()), 
		excess_0 (i_data.get_grid (1).get_excess_0 ()), 
		excess_n (i_data.get_grid (1).get_excess_n ()), 
		var_x (i_data_x),
		var_z (i_data_z),
		data_x (i_data_x.ptr (real_spectral)), 
		data_z (i_data_z.ptr (real_spectral)), 
		grid_n (i_data.get_grid (0)), 
		grid_m (i_data.get_grid (1)), 
		pos_n (&grid_n [0]), 
		messenger_ptr (i_messenger_ptr), 
		rhs (i_rhs) {
			TRACE ("Building laplace solver...");
			component_flags_x = &(i_data_x.component_flags);
			component_flags_z = &(i_data_z.component_flags);
			
			pressure = i_pressure;
			density = i_density;
			gamma = i_gamma;

			boundary_0 = i_boundary_0;
			boundary_n = i_boundary_n;
			
			// Resize the auxiliary arrays that will be needed in the solve
			ipiv.resize (ldn * (m + 2));
			
			id = messenger_ptr->get_id ();
			np = messenger_ptr->get_np ();
			
			if (id == 0) {
				x.resize (4 * (3) * (3) * ldn * np * np);
				xipiv.resize (2 * (3) * np * ldn);
			} else {
				x.resize (4 * (3) * (3) * ldn);
				xipiv.resize (1);
			}
			
			positions.resize (m + 6, 0.0);
			new_positions.resize (m + 6);
			pos_m = &positions [3];
			grad_pressure.resize (m);
			grad2_pressure.resize (m);
		
			for (int j = 0; j < m - 0; ++j) {
				pos_m [j] = grid_m [j];
			}
			
			// Calculate the overlapping region
			ntop = 0;
			nbot = 0;
			if (id != 0) {
				ntop = 1;
				messenger_ptr->send (3, &pos_m [excess_0 + 1], id - 1, 0);
			}
			if (id != np - 1) {
				nbot = 2;
				messenger_ptr->recv (3, &pos_m [m - excess_n], id + 1, 0);
				messenger_ptr->send (3, &pos_m [m - 4 - excess_n], id + 1, 1);
			} else {
				pos_m [m] = pos_m [m - 1] + (pos_m [m - 1] - pos_m [m - 2]);
			}
			if (id != 0) {
				messenger_ptr->recv (3, &pos_m [excess_0 - 3], id - 1, 1);
			} else {
				pos_m [-1] = pos_m [0] - (pos_m [1] - pos_m [0]);
			}

			for (int j = 1; j < m - 1; ++j)
			{
				grad_pressure [j] = log (pressure [j + 1] / pressure [j - 1]) / (pos_m [j + 1] - pos_m [j - 1]);
			}
			grad_pressure [0] = log (pressure [1] / pressure [0]) / (pos_m [1] - pos_m [0]);
			grad_pressure [m - 1] = log (pressure [m - 1] / pressure [m - 2]) / (pos_m [m - 1] - pos_m [m - 2]);

			grad2_pressure.resize (m);
			for (int j = 1; j < m - 1; ++j)
			{
				grad2_pressure [j] = (log (pressure [j + 1] / pressure [j]) / (pos_m [j + 1] - pos_m [j]) - log (pressure [j] / pressure [j - 1]) / (pos_m [j] - pos_m [j - 1])) / (pos_m [j + 1] - pos_m [j - 1]);
			}
			grad2_pressure [0] = grad2_pressure [1];
			grad2_pressure [m - 1] = grad2_pressure [m - 2];
			
			// Resize the main arrays for the solve
			int kl = 2;
			int ku = 1;
			int lda = 2 * kl + ku + 1;
			
			matrix.resize ((lda) * (m + 2 + kl + ku) * ldn);
			bufferl.resize (np * (m + 2) * kl * ldn);
			bufferr.resize (np * (m + 2) * ku * ldn);
			buffer.resize (np * 4 * (ku + kl) * (ku + kl));
			
			data_temp.resize ((m + 2) * ldn);

			transform = typename std::shared_ptr <plan <datatype>> (new  transforms::horizontal <datatype> (rhs, real_real));

			oodx = grid_n.get_ood ();
			oodx2 = grid_n.get_ood2 ();
			oodz = grid_m.get_ood ();
			oodz2 = grid_m.get_ood2 ();
		}

		template <class datatype>
		void pseudo_incompressible <datatype>::factorize () {
			int info;
			TRACE ("Factorizing laplace solver...");
			
			// Calculate some relevant constants
			double scalar = 4.0 * std::acos (-1.0) * std::acos (-1.0) / (pos_n [n - 1] - pos_n [0]) / (pos_n [n - 1] - pos_n [0]);
			int kl = 2;
			int ku = 1;
			int lda = 2 * kl + ku + 1;
			
			// Define some new pointers for convenience and speed
			datatype *matrix_ptr, *npos_m = &pos_m [excess_0];
			new_pos = &new_positions [3];
			
			// Update new_pos to contain the midpoints of pos_m
			for (int j = -3; j < m + 2; ++j) {
				new_pos [j] = (pos_m [j] + pos_m [j + 1]) / 2.0;
			}
			if (id == 0) {
				new_pos [-1] = 2.0 * pos_m [0] - (pos_m [0] + pos_m [1]) / 2.0;
			}
			if (id == np - 1) {
				new_pos [m - 1] = 2.0 * pos_m [m - 1] - (pos_m [m - 1] + pos_m [m - 2]) / 2.0;
				new_pos [m] = pos_m [m - 1] + (pos_m [m - 1] - new_pos [m - 3]);
			}
			new_pos += excess_0;
			
			// Generate the matrix
			#pragma omp parallel for
			for (int i = 0; i < ldn; ++i) {
				datatype *matrix_ptr = &matrix [(i) * (m + 2 + kl + ku) * lda + kl + ku + (kl + ku + excess_0) * lda];
				for (int j = 0; j < m + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0); ++j) {
					// j is the gridpoint location of the solve on the velocity grid (I think)
					matrix_ptr [(j - 2) * lda + 2] = 1.0 / (new_pos [j - 1] - new_pos [j - 2]) / (npos_m [j + 1] - npos_m [j - 1]);
					matrix_ptr [(j - 1) * lda + 1] = 
						-1.0 / (new_pos [j - 1] - new_pos [j - 2]) / (npos_m [j + 1] - npos_m [j - 1])
						- scalar * (i / 2) * (i / 2) / 2.0
						- grad_pressure [excess_0 + j] / gamma * grad_pressure [excess_0 + j] / gamma / 2.
						- grad2_pressure [excess_0 + j] / gamma / 2.;
					matrix_ptr [(j) * lda] = 
						-1.0 / (new_pos [j + 1] - new_pos [j]) / (npos_m [j + 1] - npos_m [j - 1])
						- scalar * (i / 2) * (i / 2) / 2.0
						- grad_pressure [excess_0 + j] / gamma * grad_pressure [excess_0 + j] / gamma / 2.
						- grad2_pressure [excess_0 + j] / gamma / 2.;
					matrix_ptr [(j + 1) * lda - 1] = 1.0 / (new_pos [j + 1] - new_pos [j]) / (npos_m [j + 1] - npos_m [j - 1]);
				}
				if (id == 0) {
					matrix_ptr [-3 * lda + 2] = 0.0;
					matrix_ptr [-2 * lda + 1] = 0.0;
					matrix_ptr [-1 * lda] = 1.0;
					matrix_ptr [-1] = -1.0;
					
					matrix_ptr [-2 * lda + 2] = 0.0;
					matrix_ptr [-1 * lda + 1] = 1.0 / (new_pos [0] - new_pos [-1]) / (npos_m [1] - npos_m [0]);
					matrix_ptr [0] = -1.0 / (new_pos [1] - new_pos [0]) / (npos_m [1] - npos_m [0]) - 1.0 / (new_pos [0] - new_pos [-1]) / (npos_m [1] - npos_m [0]) - scalar * (i / 2) * (i / 2);
					matrix_ptr [1 * lda - 1] = 1.0 / (new_pos [1] - new_pos [0]) / (npos_m [1] - npos_m [0]);
				}
				if (id == np - 1) {
					int j = m + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0);
					matrix_ptr [(j - 2) * lda + 2] = -1.0;
					matrix_ptr [(j - 1) * lda + 1] = 1.0;
					matrix_ptr [j * lda] = 1.0e-10;
					matrix_ptr [(j + 1) * lda - 1] = 0.0;
				}
			}
			
			linalg::block::banded_factorize (id, np, m + (nbot == 0 ? 1 : -nbot - excess_n - 1) + (id == 0 ? 1: -excess_0 - ntop), kl, ku, &matrix [(id == 0 ? 0 : 1 + excess_0) * lda], &ipiv [0], &x [0], &xipiv [0], &bufferl [0], &bufferr [0], &buffer [0], &info, ldn, lda, m + 2 + kl + ku);
		}
		
		template <class datatype>
		void pseudo_incompressible <datatype>::execute () {
			solver <datatype>::execute ();

			int info;
			TRACE ("Solving...");
			int kl = 2;
			int ku = 1;
			int lda = 2 * kl + ku + 1;
				
			linalg::scale ((m + 2) * ldn, 0.0, &data_temp [0]);

			datatype *rhs_real = rhs.ptr ();

			#pragma omp parallel for
			for (int i = 0; i < n; ++i)
			{
				int p1 = (i - 1) % n;
				int m1 = (i - 1 + n) % n;
				int g = 0;
				for (int j = 1; j < m - 1; ++j)
				{
					g = i * m + j;
					rhs_real [g] -= ((data_in [p1 * m + j] - data_in [m1 * m + j]) * oodx2 [i]) * log (density [p1 * m + j] / density [m1 * m + j]) * oodx2 [i];
					rhs_real [g] -= (data_in [g] * grad_pressure [j] / gamma + (data_in [g + 1] - data_in [g - 1]) * oodz2 [j]) * log (density [g + 1] / density [g - 1]) * oodz2 [j];
					rhs_real [g] -= data_in [g] * ((log (density [g + 1] / density [g]) * oodz [j] - log (density [g] / density [g - 1]) * oodz [j - 1]) * oodz2 [j] + (log (density [p1 * m + j] / density [g]) * oodx [i] - log (density [g] / density [m1 * m + j]) * oodx [m1]) * oodx2 [i]) * 2.;
				}
			}

			transform->execute ();

			linalg::matrix_copy (m, ldn, rhs_real, &data_temp [1], m, m + 2);
			linalg::matrix_add_scaled (m, ldn, (gamma - 1.0) / gamma, rhs.ptr (real_spectral), &data_temp [1], m, m + 2);
		
			datatype scalar = acos (-1.0) * 2.0 / (pos_n [n - 1] - pos_n [0]);
			datatype *data_ptr = &data_temp [1];
			
			datatype *npos_m = &pos_m [excess_0], *ndata_z = &data_z [excess_0], *ndata_x = &data_x [excess_0];
			
			data_ptr += excess_0;
					
			// Calculate the horizontal derivatives of the horizontal component
			#pragma omp parallel for
			for (int i = 2; i < ldn; i += 2) {
				for (int j = -1; j < m + 1 + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0); ++j) {
					data_ptr [i * (m + 2) + j] = -scalar * (i / 2) * ndata_x [(i + 1) * m + j];
					data_ptr [(i + 1) * (m + 2) + j] = scalar * (i / 2) * ndata_x [i * m + j];
				}
			}
			
			// Calculate the vertical derivatives of the vertical component
			#pragma omp parallel for
			for (int i = 0; i < ldn; ++i) {
				if (id == 0) {
					data_ptr [i * (m + 2) - 1] = 0.0;
					data_ptr [i * (m + 2)] += (ndata_z [i * m + 1] - ndata_z [i * m]) / (npos_m [1] - npos_m [0]) - ndata_z [i * m] / gamma * grad_pressure [0];
				}
				for (int j = (id == 0 ? 1 : 0); j < m + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0); ++j) {
					data_ptr [i * (m + 2) + j] += (ndata_z [i * m + j + 1] - ndata_z [i * m + j - 1]) / (npos_m [j + 1] - npos_m [j - 1]) - ndata_z [i * m + j] / gamma * grad_pressure [j];
				}
				if (id == np - 1) {
					data_ptr [i * (m + 2) + (m + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0))] = 0.0;
				}
			}
		
			linalg::block::banded_solve (id, np, m + (nbot == 0 ? 1 : -nbot - excess_n - 1) + (id == 0 ? 1: -excess_0 - ntop), kl, ku, &matrix [(id == 0 ? 0 : 1 + excess_0) * 6], &ipiv [0], &data_temp [(id == 0 ? 0 : 1 + excess_0)], &x [0], &xipiv [0], &bufferl [0], &bufferr [0], &info, ldn, lda, m + 2 + kl + ku, m + 2);
			
			// We can define our frame to be one with no net vertical flux
			linalg::scale (2 * (m + 2), 0.0, &data_temp [0]);
			
			// Communicate the edges of the pressure
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

			ndata_z = &data_z [excess_0];
			ndata_x = &data_x [excess_0];
			
			// Update the velocities with the pressure derivatives
			#pragma omp parallel for
			for (int i = 2; i < ldn; i += 2) {
				for (int j = 0; j < m + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0); ++j) {
					ndata_x [i * m + j] += scalar * (i / 2) * (data_ptr [(i + 1) * (m + 2) + j] + data_ptr [(i + 1) * (m + 2) + j - 1]) / 2.0;
					ndata_x [(i + 1) * m + j] -= scalar * (i / 2) * (data_ptr [i * (m + 2) + j] + data_ptr [i * (m + 2) + j - 1]) / 2.0;
				}
			}

			#pragma omp parallel for
			for (int i = 2; i < ldn; ++i) {
				for (int j = 0; j < m + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0); ++j) {
					ndata_z [i * m + j] -= (data_ptr [i * (m + 2) + j] - data_ptr [i * (m + 2) + j - 1]) / (new_pos [j] - new_pos [j - 1])
						- (data_ptr [i * (m + 2) + j] + data_ptr [i * (m + 2) + j - 1]) / 2.0 / gamma * grad_pressure [j];
				}
				for (int j = 0; j < m + 1 - (nbot == 0 ? 0 : excess_n + 1) - excess_0; ++j) {
					data_out [i * m + j + excess_0] = (data_temp [i * (m + 2) + j + excess_0] + data_temp [i * (m + 2) + j + 1 + excess_0]) / 2.0;
				}
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

#ifdef CHECKNAN
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < ldn; ++i) {
					if (std::isnan (data_out [i * m + j])) {
						FATAL ("Nan in laplace solver.");
						throw linalg::exceptions::nan ();
					}
				}
			}
#endif

			var_x.component_flags &= ~grids::updated;
			var_x.last_update = get_state ();
			var_z.component_flags &= ~grids::updated;
			var_z.last_update = get_state ();

			TRACE ("Solved");
		}

		template class pseudo_incompressible <double>;
	} /* solvers */
} /* plans */
