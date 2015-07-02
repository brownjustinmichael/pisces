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
		pseudo_incompressible <datatype>::pseudo_incompressible (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, mpi::messenger* i_messenger_ptr, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_0, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_n, datatype* i_data, datatype *i_data_x, datatype *i_data_z, int *i_element_flags, int *i_component_flags, int * i_component_flags_x, int *i_component_flags_z) : solver <datatype> (i_element_flags, i_component_flags), n (i_grid_n.get_n ()), ldn (i_grid_n.get_ld ()), m (i_grid_m.get_n ()), excess_0 (i_grid_m.get_excess_0 ()), excess_n (i_grid_m.get_excess_n ()), data (i_data), data_x (i_data_x), data_z (i_data_z), grid_n (i_grid_n), grid_m (i_grid_m), pos_n (&grid_n [0]), messenger_ptr (i_messenger_ptr) {
			TRACE ("Building laplace solver...");
			component_flags_x = i_component_flags_x;
			component_flags_z = i_component_flags_z;
			
			boundary_0 = i_boundary_0;
			boundary_n = i_boundary_n;
			
			// Resize the auxiliary arrays that will be needed in the solve
			ipiv.resize (ldn * (m + 2));
			
			id = messenger_ptr->get_id ();
			np = messenger_ptr->get_np ();
			
			// if (id == 0) {
				x.resize (2 * m * ldn + 8 * np * ldn, 0.0);
				xipiv.resize (2 * np * ldn);
			// } else {
				// x.resize (4 * (3) * (3) * ldn);
				// xipiv.resize (1);
			// }
			
			positions.resize (m + 2);
			pos_m = &positions [1];
			diff.resize (m + 1);
			diff2.resize (m);
		
			for (int j = 0; j < m - 0; ++j) {
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

			matrix.resize (4 * m * ldn);
		}

		template <class datatype>
		void pseudo_incompressible <datatype>::factorize () {
			int info;
			TRACE ("Factorizing laplace solver...");
			
			// Calculate some relevant constants
			double scalar = 4.0 * std::acos (-1.0) * std::acos (-1.0) / (pos_n [n - 1] - pos_n [0]) / (pos_n [n - 1] - pos_n [0]);

			// Define some new pointers for convenience and speed
			datatype *sub_ptr = &matrix [0], *diag_ptr = &matrix [m * ldn], *sup_ptr = &matrix [2 * m * ldn];
			
			// Generate the matrix
			for (int i = 0; i < ldn; ++i) {
				for (int j = 0; j < m; ++j) {
					sup_ptr [j] = 1.0 / diff [j + 1] / diff2 [j];
					diag_ptr [j] = (-1.0 / diff [j + 1] - 1.0 / diff [j]) / diff2 [j] - scalar * (i / 2) * (i / 2);
					sub_ptr [j] = 1.0 / diff [j] / diff2 [j];
				}
				if (id == 0) {
					sup_ptr [0] = 2.0 / diff [1] / diff2 [0];
					diag_ptr [0] = -2.0 / diff [1] / diff2 [0] - scalar * (i / 2) * (i / 2);
				}
				if (id == np - 1) {
					DEBUG (diff2 [m - 1] << diff [m - 1]);
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
			int info;
			static int count = 0;
			TRACE ("Solving...");
			bool retransform = false;
			// No net vertical flux
			// linalg::scale (2 * m, 0.0, data_z);
			++count;
			if (count % 2 == 0) return;

			// if (!(*component_flags_x & transformed_vertical)) {
				linalg::scale (m * ldn, 0.0, &data [0]);

				datatype scalar = acos (-1.0) * 2.0 / (pos_n [n - 1] - pos_n [0]);
				
				// Calculate the horizontal derivatives of the horizontal component
				for (int i = 2; i < ldn; i += 2) {
					for (int j = 0; j < m; ++j) {
						data [i * m + j] += -scalar * (i / 2) * data_x [(i + 1) * m + j];
						data [(i + 1) * m + j] += scalar * (i / 2) * data_x [i * m + j];
					}
				}
				
				// Calculate the vertical derivatives of the vertical component
				for (int i = 0; i < ldn; ++i) {
					if (id == 0) {
						data [i * m] += (data_z [i * m + 1] - data_z [i * m]) / diff [1];
					}
					for (int j = 1; j < m - 1; ++j) {
						data [i * m + j] += (data_z [i * m + j + 1] - data_z [i * m + j - 1]) / diff2 [j];
					}
					if (id == np - 1) {
						data [i * m + m - 1] = (data_z [i * m + m - 1] - data_z [i * m + m - 2]) / diff [m - 1];
					}
				}

				std::stringstream debug;
				for (int j = 0; j < 5; ++j)
				{
					DEBUG (pos_m [0] << " " << pos_m [1] << pos_m [2])
					for (int i = 0; i < ldn; ++i)
					{
						debug << data [i * m + j] << " ";
					}
					DEBUG (debug.str ());
					debug.str ("");
				}

				for (int j = m - 5; j < m; ++j)
				{
					DEBUG (pos_m [m - 1] << " " << pos_m [m - 2] << " " << pos_m [m - 3])
					for (int i = 0; i < ldn; ++i)
					{
						debug << data [i * m + j] << " ";
					}
					DEBUG (debug.str ());
					debug.str ("");
				}

				linalg::block::tridiag_solve (id, np, m - (excess_0 ? excess_0 + 1 : 0) - (excess_n ? excess_n + 2 : 0), &matrix [excess_0], &matrix [excess_0 + m * ldn], &matrix [excess_0 + 2 * m * ldn], &matrix [excess_0 + 3 * m * ldn], &ipiv [0], &data [excess_0], &x [0], &xipiv [0], &info, ldn, m, m);
				
				// We can define our frame to be one with no net vertical flux
				linalg::scale (2 * m, 0.0, &data [0]);

				// Communicate the edges of the pressure
				if (boundary_n) {
					boundary_n->send (&data [m - 2 * excess_n - 1], m, excess_n);
				}
				if (boundary_0) {
					boundary_0->receive (&data [0], m, excess_0, 0.0);
					boundary_0->send (&data [excess_0], m, excess_0 + 1);
				}
				if (boundary_n) {
					boundary_n->receive (&data [m - excess_n - 1], m, excess_n + 1, 0.0);
				}
				
				// Update the velocities with the pressure derivatives
				for (int i = 2; i < ldn; ++i) {
					for (int j = 1; j < m - 1; ++j) {
						data_z [i * m + j] -= (data [i * m + j + 1] - data [i * m + j - 1]) / diff2 [j];
					}
				}
				
				for (int i = 2; i < ldn; i += 2) {
					for (int j = 0; j < m; ++j) {
						data_x [i * m + j] += scalar * (i / 2) * data [(i + 1) * m + j];
						data_x [(i + 1) * m + j] -= scalar * (i / 2) * data [i * m + j];
					}
				}

				if (boundary_n) {
					boundary_n->send (&data_z [m - 2 * excess_n - 1], m, excess_n);
				}
				if (boundary_0) {
					boundary_0->receive (&data_z [0], m, excess_0, 0.0);
					boundary_0->send (&data_z [excess_0], m, excess_0 + 1);
				}
				if (boundary_n) {
					boundary_n->receive (&data_z [m - excess_n - 1], m, excess_n + 1, 0.0);
				}

				if (boundary_n) {
					boundary_n->send (&data_x [m - 2 * excess_n - 1], m, excess_n);
				}
				if (boundary_0) {
					boundary_0->receive (&data_x [0], m, excess_0, 0.0);
					boundary_0->send (&data_x [excess_0], m, excess_0 + 1);
				}
				if (boundary_n) {
					boundary_n->receive (&data_x [m - excess_n - 1], m, excess_n + 1, 0.0);
				}

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

			TRACE ("Solved");
		}

		template class pseudo_incompressible <double>;
	} /* solvers */
} /* plans */
