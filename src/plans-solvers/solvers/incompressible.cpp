/*!**********************************************************************
 * \file solver_two_d.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-14.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cmath>
#include <sstream>

#include "linalg/utils.hpp"
#include "linalg/linalg.hpp"
#include "linalg/interpolate.hpp"
#include "linalg/exceptions.hpp"
#include "linalg-block/banded.hpp"

#include "incompressible.hpp"
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
	incompressible_corrector <datatype>::incompressible_corrector (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, mpi::messenger* i_messenger_ptr, datatype *i_rhs, datatype* i_data, datatype *i_data_x, datatype *i_data_z, int *i_element_flags, int *i_component_flags, int * i_component_flags_x, int *i_component_flags_z) :
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
		pos_m = &positions [3];
		
		for (int j = 0; j < m - 0; ++j) {
			pos_m [j] = grid_m [j];
		}
		
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
		
		int kl = 2;
		int ku = 1;
		int lda = 2 * kl + ku + 1;
		
		matrix.resize ((lda) * (m + 2 + kl + ku) * ldn);
		bufferl.resize (np * (m + 2) * kl * ldn);
		bufferr.resize (np * (m + 2) * ku * ldn);

		data_temp.resize ((m + 2) * ldn);
	}

	template <class datatype>
	void incompressible_corrector <datatype>::factorize () {
		std::stringstream debug;
		int info;
		TRACE ("Factorizing laplace solver...");

		double scalar = 4.0 * std::acos (-1.0) * std::acos (-1.0) / (pos_n [n - 1] - pos_n [0]) / (pos_n [n - 1] - pos_n [0]);
		std::vector <datatype> positions (m + 2 + 6);
		datatype *matrix_ptr, *new_pos = &positions [3], *npos_m = &pos_m [excess_0];
		int kl = 2;
		int ku = 1;
		int lda = 2 * kl + ku + 1;
		
		for (int j = -3; j < m + 3; ++j) {
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
		
		for (int i = 0; i < ldn; ++i) {
			matrix_ptr = &matrix [(i) * (m + 2 + kl + ku) * lda + kl + ku + (kl + ku + excess_0) * lda];
			for (int j = 0; j < m + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0); ++j) {
				// j is the gridpoint location of the solve on the velocity grid (I think)
				matrix_ptr [(j - 2) * lda + 2] = 1.0 / (new_pos [j - 1] - new_pos [j - 2]) / (npos_m [j + 1] - npos_m [j - 1]);
				matrix_ptr [(j - 1) * lda + 1] = -1.0 / (new_pos [j - 1] - new_pos [j - 2]) / (npos_m [j + 1] - npos_m [j - 1]) - scalar * (i / 2) * (i / 2) / 2.0;
				matrix_ptr [(j) * lda] = -1.0 / (new_pos [j + 1] - new_pos [j]) / (npos_m [j + 1] - npos_m [j - 1]) - scalar * (i / 2) * (i / 2) / 2.0;
				matrix_ptr [(j + 1) * lda - 1] = 1.0 / (new_pos [j + 1] - new_pos [j]) / (npos_m [j + 1] - npos_m [j - 1]);
			}
			if (id == 0) {
				matrix_ptr [-3 * lda + 2] = 0.0;
				matrix_ptr [-2 * lda + 1] = 0.0;
				matrix_ptr [-1 * lda] = 1.0;
				matrix_ptr [-1] = -1.0;
				// matrix_ptr [-15] = 0.0;
				// matrix_ptr [10] = 1.0 / (new_pos [-1] - new_pos [-2]) / (npos_m [0] - npos_m [-1]);
				// matrix_ptr [-5] = -1.0 / (new_pos [0] - new_pos [-1]) / (npos_m [0] - npos_m [-1]) - 1.0 / (new_pos [-1] - new_pos [-2]) / (npos_m [0] - npos_m [-1]) - scalar * (i / 2) * (i / 2);
				// matrix_ptr [0] = 1.0 / (new_pos [0] - new_pos [-1]) / (npos_m [0] - npos_m [-1]);

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
				matrix_ptr [j * lda - 1] = 0.0;
				// matrix_ptr [(j - 2) * 6 + 3] = 1.0 / (new_pos [j - 1] - new_pos [j - 2]) / (npos_m [j + 1] - npos_m [j - 1]);
				// matrix_ptr [(j - 1) * 6 + 2] = -1.0 / (new_pos [j - 1] - new_pos [j - 2]) / (npos_m [j + 1] - npos_m [j - 1]) - scalar * (i / 2) * (i / 2) / 2.0;
				// matrix_ptr [(j) * 6 + 1] = -1.0 / (new_pos [j + 1] - new_pos [j]) / (npos_m [j + 1] - npos_m [j - 1]) - scalar * (i / 2) * (i / 2) / 2.0;
				// matrix_ptr [(j + 1) * 6] = 1.0 / (new_pos [j + 1] - new_pos [j]) / (npos_m [j + 1] - npos_m [j - 1]);
			}
		}
		
		// throw 0;
		linalg::p_block_banded_factorize (id, np, m + (nbot == 0 ? 1 : -nbot - excess_n - 1) + (id == 0 ? 1: -excess_0 - ntop), kl, ku, &matrix [(id == 0 ? 0 : 1 + excess_0) * lda], &ipiv [0], &x [0], &xipiv [0], &bufferl [0], &bufferr [0], &info, ldn, lda, m + 2 + kl + ku);
		// throw 0;
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
		
		linalg::scale ((m + 2) * ldn, 0.0, &data_temp [0]);
		
		std::shared_ptr <plans::plan <datatype> > transform_x = std::shared_ptr <plans::plan <datatype> > (new plans::vertical_transform <datatype> (n, m, data_x, NULL, 0x00, element_flags, component_flags_x));
		std::shared_ptr <plans::plan <datatype> > transform_z = std::shared_ptr <plans::plan <datatype> > (new plans::vertical_transform <datatype> (n, m, data_z, NULL, 0x00, element_flags, component_flags_z));
		
		if (*component_flags_x & transformed_vertical) {
			transform_x->execute ();
			transform_z->execute ();
			retransform = true;
		}

		if (!(*component_flags_x & transformed_vertical)) {
			DEBUG ("Here");
			datatype scalar = acos (-1.0) * 2.0 / (pos_n [n - 1] - pos_n [0]);
			datatype *data_ptr = &data_temp [1];

			std::vector <datatype> positions (m + 2 + 6);
			datatype *new_pos = &positions [3], *npos_m = &pos_m [excess_0], *ndata_z = &data_z [excess_0], *ndata_x = &data_x [excess_0];

			for (int j = -3; j < m + 3; ++j) {
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
			data_ptr += excess_0;
		
			
			for (int i = 2; i < ldn; i += 2) {
				for (int j = -1; j < m + 1 + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0); ++j) {
					data_ptr [i * (m + 2) + j] = -scalar * (i / 2) * ndata_x [(i + 1) * m + j];
					data_ptr [(i + 1) * (m + 2) + j] = scalar * (i / 2) * ndata_x [i * m + j];
				}
			}
			
			for (int i = 0; i < ldn; ++i) {
				if (id == 0) {
					data_ptr [i * (m + 2) - 1] = 0.0;
					data_ptr [i * (m + 2)] += (ndata_z [i * m + 1] - ndata_z [i * m]) / (npos_m [1] - npos_m [0]);
				}
				for (int j = (id == 0 ? 1 : 0); j < m + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0); ++j) {
					data_ptr [i * (m + 2) + j] += (ndata_z [i * m + j + 1] - ndata_z [i * m + j - 1]) / (npos_m [j + 1] - npos_m [j - 1]);
				}
				if (id == np - 1) {
					data_ptr [i * (m + 2) + (m + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0))] = 0.0;
				}
			}
			linalg::p_block_banded_solve (id, np, m + (nbot == 0 ? 1 : -nbot - excess_n - 1) + (id == 0 ? 1: -excess_0 - ntop), kl, ku, &matrix [(id == 0 ? 0 : 1 + excess_0) * 6], &ipiv [0], &data_temp [(id == 0 ? 0 : 1 + excess_0)], &x [0], &xipiv [0], &bufferl [0], &bufferr [0], &info, ldn, lda, m + 2 + kl + ku, m + 2);
			// throw 0;

			linalg::scale (2 * (m + 2), 0.0, &data_temp [0]);

			for (int i = 2; i < ldn; ++i) {
				for (int j = 0; j < m + 1 - (nbot == 0 ? 0 : excess_n + 1) - excess_0; ++j) {
					data [i * m + j + excess_0] = (data_temp [i * (m + 2) + j + excess_0] + data_temp [i * (m + 2) + j + 1 + excess_0]) / 2.0;
				}
			}

			for (int i = 2; i < ldn; ++i) {
				for (int j = 1; j < m - 1 + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0); ++j) {
					ndata_z [i * m + j] -= (data_ptr [i * (m + 2) + j] - data_ptr [i * (m + 2) + j - 1]) / (new_pos [j] - new_pos [j - 1]);
				}
			}

			for (int i = 2; i < ldn; i += 2) {
				for (int j = 0; j < m + (nbot == 0 ? 0 : -excess_n - 1) + (id == 0 ? 0: -excess_0); ++j) {
					ndata_x [i * m + j] += scalar * (i / 2) * (data_ptr [(i + 1) * (m + 2) + j] + data_ptr [(i + 1) * (m + 2) + j - 1]) / 2.0;
					ndata_x [(i + 1) * m + j] -= scalar * (i / 2) * (data_ptr [i * (m + 2) + j] + data_ptr [i * (m + 2) + j - 1]) / 2.0;
				}
			}
			
			for (int i = 0; i < n; ++i) {
				for (int j = excess_0 - 1; j >= 0; --j) {
					data_x [i * m + j] = data_x [i * m + j + 1] + (data_x [i * m + j + 2] - data_x [i * m + j + 1]) / (pos_m [j + 2] - pos_m [j + 1]) * (pos_m [j] - pos_m [j + 1]);
					data_z [i * m + j] = data_z [i * m + j + 1] + (data_z [i * m + j + 2] - data_z [i * m + j + 1]) / (pos_m [j + 2] - pos_m [j + 1]) * (pos_m [j] - pos_m [j + 1]);
				}
				
				for (int j = m - excess_n - 1; j < m; ++j) {
					data_x [i * m + j] = data_x [i * m + j - 1] + (data_x [i * m + j - 2] - data_x [i * m + j - 1]) / (pos_m [j - 2] - pos_m [j - 1]) * (pos_m [j] - pos_m [j - 1]);
					data_z [i * m + j] = data_z [i * m + j - 1] + (data_z [i * m + j - 2] - data_z [i * m + j + 1]) / (pos_m [j + 2] - pos_m [j - 1]) * (pos_m [j] - pos_m [j - 1]);
				}
			}
			
			// for (int i = 0; i < n; ++i) {
			// 	data [i * m] = 0.0;
			// 	data [i * m + m - 1] = 0.0;
			// }
			// for (int j = 1; j < m - 1; ++j) {
			// 	for (int i = 0; i < n; i += 2) {
			// 		data [i * m + j] = -scalar * (i / 2) * data_x [(i + 1) * m + j] + (data_z [i * m + j + 1] - data_z [i * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]);
			// 		data [(i + 1) * m + j] = scalar * (i / 2) * data_x [i * m + j] + (data_z [(i + 1) * m + j + 1] - data_z [(i + 1) * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]);
			// 	}
			// }

			linalg::scale (2 * m, 0.0, data_z);
			linalg::scale (2 * m, 0.0, data_x);

			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < ldn; ++i) {
					if (std::isnan (data [i * m + j])) {
						FATAL ("Nan in laplace solver.");
						throw linalg::exceptions::nan ();
					}
				}
			}
		} else {
			FATAL ("SHOULDN'T BE HERE");
			throw 0;
		}
		
		if (retransform) {
			transform_x->execute ();
			transform_z->execute ();
		}

		TRACE ("Solved");
	}

	template class incompressible_corrector <double>;
} /* plans */
