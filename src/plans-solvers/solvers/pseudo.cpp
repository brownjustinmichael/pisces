/*!**********************************************************************
 * \file pseudo_incompressible.cpp
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
#include "linalg-block/tridiagonal.hpp"

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
		pseudo_incompressible <datatype>::pseudo_incompressible (mpi::messenger* i_messenger_ptr, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_0, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_n, datatype i_gamma, grids::variable <datatype> &i_data, grids::variable <datatype> &i_data_out, grids::variable <datatype> &i_data_x, grids::variable <datatype> &i_data_z, datatype *i_pressure) : 
		solver <datatype> (i_data, i_data_out, this->get_state_in (), this->get_state ()), 
		n (i_data.get_grid (0).get_n ()), 
		ldn (i_data.get_grid (0).get_ld ()), 
		m (i_data.get_grid (1).get_n ()), 
		excess_0 (i_data.get_grid (1).get_excess_0 ()), 
		excess_n (i_data.get_grid (1).get_excess_n ()), 
		gamma (i_gamma), 
		var_x (i_data_x),
		var_z (i_data_z),
		data_x (i_data_x.ptr ()), 
		data_z (i_data_z.ptr ()), 
		grid_n (i_data.get_grid (0)), 
		grid_m (i_data.get_grid (1)), 
		pos_n (&grid_n [0]), 
		messenger_ptr (i_messenger_ptr) {
			TRACE ("Building laplace solver...");
			component_flags_x = &(i_data_x.component_flags);
			component_flags_z = &(i_data_z.component_flags);

			pressure = i_pressure;
			
			boundary_0 = i_boundary_0;
			boundary_n = i_boundary_n;
			
			// Resize the auxiliary arrays that will be needed in the solve
			ipiv.resize (ldn * (m + 2));
			
			id = messenger_ptr->get_id ();
			np = messenger_ptr->get_np ();
			
			x.resize (2 * m * ldn + 8 * np * ldn, 0.0);
			xipiv.resize (2 * np * ldn);
			
			positions.resize (m + 2);
			pos_m = &positions [1];
			diff.resize (m + 1);
			diff2.resize (m);
			grad_pressure.resize (m);
		
			for (int j = 0; j < m; ++j) {
				pos_m [j] = grid_m [j];
			}

			// Calculate the overlapping region
			if (id != 0) {
				ntop = 1;
				messenger_ptr->send (1, &pos_m [excess_0 + 1], id - 1, 0);
			}
			if (id != np - 1) {
				nbot = 2;
				messenger_ptr->recv (1, &pos_m [m - excess_n], id + 1, 0);
				messenger_ptr->send (1, &pos_m [m - 2 - excess_n], id + 1, 1);
			} else {
				pos_m [m] = pos_m [m - 1] + (pos_m [m - 1] - pos_m [m - 2]);
			}
			if (id != 0) {
				messenger_ptr->recv (1, &pos_m [excess_0 - 1], id - 1, 1);
			} else {
				pos_m [-1] = pos_m [0] - (pos_m [1] - pos_m [0]);
			}
			
			for (int j = 0; j < m; ++j)
			{
				diff [j] = pos_m [j] - pos_m [j - 1];
				diff2 [j] = pos_m [j + 1] - pos_m [j - 1];
			}
			diff [m] = pos_m [m] - pos_m [m - 1];

			for (int j = 1; j < m - 1; ++j)
			{
				grad_pressure [j] = (pressure [j + 1] - pressure [j - 1]) / diff2 [j];
			}
			grad_pressure [0] = (pressure [1] - pressure [0]) / diff [1];
			grad_pressure [m - 1] = (pressure [m - 1] - pressure [m - 2]) / diff [m - 1];

			grad2_pressure.resize (m);
			for (int j = 1; j < m - 1; ++j)
			{
				grad2_pressure [j] = ((pressure [j + 1] - pressure [j]) / diff [j + 1] - (pressure [j] - pressure [j - 1]) / diff [j]) / diff2 [j];
			}
			grad2_pressure [0] = grad2_pressure [1];
			grad2_pressure [m - 1] = grad2_pressure [m - 2];

			matrix.resize (4 * m * ldn);
		}

		template <class datatype>
		void pseudo_incompressible <datatype>::factorize () {
			int info;
			TRACE ("Factorizing psueudo solver...");
			
			// Calculate some relevant constants
			double scalar = 4.0 * std::acos (-1.0) * std::acos (-1.0) / (pos_n [n - 1] - pos_n [0]) / (pos_n [n - 1] - pos_n [0]);

			// Define some new pointers for convenience and speed
			datatype *sub_ptr = &matrix [0], *diag_ptr = &matrix [m * ldn], *sup_ptr = &matrix [2 * m * ldn];
			
			// Generate the matrix
			for (int i = 0; i < ldn; ++i) {
				for (int j = 0; j < m; ++j) {
					sup_ptr [j] = 1.0 / diff [j + 1] / diff2 [j] + grad_pressure [j] / gamma / pressure [j] / diff2 [j];
					diag_ptr [j] = (-1.0 / diff [j + 1] - 1.0 / diff [j]) / diff2 [j] - scalar * (i / 2) * (i / 2) + scalar * (grad2_pressure [j] / gamma / pressure [j] + (1. / gamma - 2.) / pressure [j] / pressure [j] * grad_pressure [j] * grad_pressure [j]);
					sub_ptr [j] = 1.0 / diff [j] / diff2 [j] - grad_pressure [j] / gamma / pressure [j] / diff2 [j];
				}
				if (id == 0) {
					sup_ptr [0] = 2.0 / diff [1] / diff2 [0];
					diag_ptr [0] = -2.0 / diff [1] / diff2 [0] - scalar * (i / 2) * (i / 2);
				}
				if (id == np - 1) {
					diag_ptr [m - 1] = -2.0 / diff [m - 1] / diff2 [m - 1] - scalar * (i / 2) * (i / 2);
					sub_ptr [m - 1] = 2.0 / diff [m - 1] / diff2 [m - 1];
				}
				sub_ptr += m;
				diag_ptr += m;
				sup_ptr += m;
			}
			linalg::block::tridiag_factorize (id, np, m - (excess_0 ? excess_0 + 1 : 0) - (excess_n ? excess_n + 2 : 0), &matrix [excess_0], &matrix [excess_0 + m * ldn], &matrix [excess_0 + 2 * m * ldn], &matrix [excess_0 + 3 * m * ldn], &ipiv [0], &x [0], &xipiv [0], &info, ldn, m);
		}
		
		template <class datatype>
		void pseudo_incompressible <datatype>::execute () {
			solver <datatype>::execute ();

			int info;
			TRACE ("Solving...");
			bool retransform = false;
			// No net vertical flux
			static int count = 0;
			count++;
			if (count % 2) return;

			linalg::scale (m * ldn, 0.0, data_out);

			datatype scalar = acos (-1.0) * 2.0 / (pos_n [n - 1] - pos_n [0]);
			
			// Calculate the horizontal derivatives of the horizontal component
			for (int i = 2; i < ldn; i += 2) {
				for (int j = 0; j < m; ++j) {
					data_out [i * m + j] += -pressure [j] * scalar * (i / 2) * data_x [(i + 1) * m + j];
					data_out [(i + 1) * m + j] += pressure [j] * scalar * (i / 2) * data_x [i * m + j];
				}
			}
			
			// Calculate the vertical derivatives of the vertical component
			for (int i = 0; i < ldn; ++i) {
				if (id == 0) {
					data_out [i * m] += pressure [0] * (data_z [i * m + 1] - data_z [i * m]) / diff [1];
				}
				for (int j = 1; j < m - 1; ++j) {
					data_out [i * m + j] += pressure [j] * (data_z [i * m + j + 1] - data_z [i * m + j - 1]) / diff2 [j] + grad_pressure [j] * data_z [i * m + j] / gamma;
				}
				if (id == np - 1) {
					data_out [i * m + m - 1] = pressure [m - 1] * (data_z [i * m + m - 1] - data_z [i * m + m - 2]) / diff [m - 1];
				}
			}

			linalg::block::tridiag_solve (id, np, m - (excess_0 ? excess_0 + 1 : 0) - (excess_n ? excess_n + 2 : 0), &matrix [excess_0], &matrix [excess_0 + m * ldn], &matrix [excess_0 + 2 * m * ldn], &matrix [excess_0 + 3 * m * ldn], &ipiv [0], &data [excess_0], &x [0], &xipiv [0], &info, ldn, m, m);
			
			// We can define our frame to be one with no net vertical flux
			linalg::scale (2 * m, 0.0, data_out);

			// Communicate the edges of the pressure
			boundaries::boundary_match (m, boundary_0, boundary_n, data_out, excess_0, excess_n);
			
			// Update the velocities with the pressure derivatives
			for (int i = 2; i < ldn; ++i) {
				for (int j = 1; j < m - 1; ++j) {
					data_z [i * m + j] -= (data_out [i * m + j + 1] - data_out [i * m + j - 1]) / diff2 [j] / pressure [j] + data_out [j] * grad_pressure [j] / pressure [j] / pressure [j] / gamma;
				}
			}
			
			for (int i = 2; i < ldn; i += 2) {
				for (int j = 0; j < m; ++j) {
					data_x [i * m + j] += scalar * (i / 2) * data_out [(i + 1) * m + j];
					data_x [(i + 1) * m + j] -= scalar * (i / 2) * data_out [i * m + j];
				}
			}

			boundaries::boundary_match (m, boundary_0, boundary_n, data_x, excess_0, excess_n);
			boundaries::boundary_match (m, boundary_0, boundary_n, data_z, excess_0, excess_n);

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

			var_x.component_flags &= ~transforms::updated;
			var_x.last_update = get_state ();
			var_z.component_flags &= ~transforms::updated;
			var_z.last_update = get_state ();

			TRACE ("Solved");
		}

		template class pseudo_incompressible <double>;
	} /* solvers */
} /* plans */
