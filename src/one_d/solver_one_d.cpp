/*!***********************************************************************
 * \file solver_one_d.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-13.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "solver_one_d.hpp"

#include "../utils/messenger.hpp"
#include <cmath>
#include "../config.hpp"
#include "../utils/utils.hpp"
#include "../utils/solver_utils.hpp"
#include "../utils/block_solver.hpp"
#include "../utils/interpolate.hpp"
#include "element_one_d.hpp"

namespace one_d
{
	template <class datatype>
	solver <datatype>::solver (bases::grid <datatype> &i_grid, utils::messenger* i_messenger_ptr, datatype& i_timestep, datatype& i_alpha_0, datatype& i_alpha_n, datatype* i_data, int *i_element_flags, int *i_component_flags) :
	bases::solver <datatype> (i_element_flags, i_component_flags), 
	n (i_grid.get_n ()),
	ld (i_grid.get_ld ()),
	grid (i_grid),
	data (i_data),
	messenger_ptr (i_messenger_ptr),
	timestep (i_timestep), 
	alpha_0 (i_alpha_0), 
	alpha_n (i_alpha_n), 
	positions (&(grid [0])),
	excess_0 (grid.get_excess_0 ()), 
	excess_n (grid.get_excess_n ()),
	default_matrix (i_grid.get_data (0)) {
		matrix.resize (n * n);
		values_0.resize (1);
		values_n.resize (1);
		implicit_rhs_vec.resize (1 * n);
		explicit_rhs_vec.resize (1 * n);
		real_rhs_vec.resize (1 * n);
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
			messenger_ptr->template send <datatype> (excess_n, &(positions [n - excess_n]), messenger_ptr->get_id () + 1, 0);
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
		factorized_matrix.resize ((n + ex_excess_0 + ex_excess_n + nbot + ntop) * (n + ex_excess_0 + ex_excess_n + nbot + ntop), 0.0);
		ipiv.resize (n); // TODO Should be n - ntop - nbot - excess_0 - excess_n
		data_temp.resize ((n + ex_excess_0 + ex_excess_n + nbot + ntop) * (1));
	}

	template <class datatype>
	void solver <datatype>::_factorize () {
		int info, lda = n + ex_excess_0 + ex_excess_n + nbot + ntop;
		TRACE ("Factorizing..." << messenger_ptr->get_id ());

		utils::matrix_scale (lda, lda, 0.0, &factorized_matrix [0], lda);
		utils::matrix_copy (n, n, default_matrix, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1)], n, lda);

		utils::matrix_add_scaled (n - excess_n - excess_0 - 2, n, timestep, &matrix [0] + excess_0 + 1, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + 1 + excess_0], n, lda);
		if (ntop != 0) {
			utils::matrix_add_scaled (ntop, n, alpha_0 * timestep, &matrix [0] + excess_0, &factorized_matrix [(ntop + ex_excess_0) * lda], n, n + ex_excess_0 + ex_excess_n + ntop + nbot);
			utils::interpolate (ex_excess_0, n, n, timestep, 1.0, positions, &matrix [0], &positions_0 [0], &factorized_matrix [(ntop + ex_excess_0) * lda + ntop], n, lda);
			utils::matrix_add_scaled (ntop, n, alpha_0 * timestep, &matrix [0] + excess_0, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + excess_0], n, n + ex_excess_0 + ex_excess_n + ntop + nbot);
		}
		if (nbot != 0) {
			utils::matrix_add_scaled (nbot, n, alpha_n * timestep, &matrix [0] + n - nbot - excess_n, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + n - nbot - excess_n], n, lda);
			utils::interpolate (ex_excess_n, n, n, timestep, 1.0, positions, &matrix [0], &positions_n [0], &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + n], n, lda);
			utils::matrix_add_scaled (nbot, n, alpha_n * timestep, &matrix [0] + n - nbot - excess_n, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + n + ex_excess_n], n, lda);
		}
		
		utils::p_block_matrix_factorize (messenger_ptr->get_id (), messenger_ptr->get_np (), n - excess_0 - excess_n - ntop - nbot, excess_0 + ex_excess_0 + 2 * ntop, excess_n + ex_excess_n + 2 * nbot, &factorized_matrix [0], &ipiv [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, lda, sqrt ((int) boundary_matrix.size ()));

		TRACE ("Done.");
	}

	template <class datatype>
	void solver <datatype>::_solve () {
		int info, lda = n + ex_excess_0 + ex_excess_n + nbot + ntop;
		TRACE ("Executing solve...");
		utils::scale (lda, 0.0, &data_temp [0]);
	
		utils::scale ((1) * lda, 0.0, &data_temp [0]);
				
		if (!(*component_flags & first_run)) {
			utils::copy (1, &data [0], &values_0 [0], n);
			utils::copy (1, &data [n - 1], &values_n [0], n);
			*component_flags |= first_run;
		}
		
		utils::matrix_add_scaled (n - excess_0 - excess_n, 1, timestep, &implicit_rhs_vec [excess_0], &data_temp [ex_excess_0 + ntop + excess_0], n, lda);
		utils::matrix_add_scaled (n - excess_0 - excess_n, 1, timestep, &explicit_rhs_vec [excess_0], &data_temp [ex_excess_0 + ntop + excess_0], n, lda);
						
		if (ntop != 0) {
			utils::interpolate (ex_excess_0, 1, n - excess_0 - excess_n, 1.0, 1.0, positions + excess_0, &data_temp [ex_excess_0 + ntop + excess_0], &positions_0 [0], &data_temp [ntop], lda, lda);
			utils::scale (1, alpha_0, &data_temp [0] + ntop + ex_excess_0 + excess_0, lda);
			utils::copy (1, &data_temp [ntop + ex_excess_0 + excess_0], &data_temp [0], lda, lda);
		} else {
			utils::copy (1, &values_0 [0], &data_temp [0], 1, lda);
		}
		if (nbot != 0) {
			utils::interpolate (ex_excess_n, 1, n - excess_0 - excess_n, 1.0, 1.0, positions + excess_0, &data_temp [ex_excess_0 + ntop + excess_0], &positions_n [0], &data_temp [lda - nbot - ex_excess_n], lda, lda);
			utils::scale (1, alpha_n, &data_temp [0] + lda - 2 * nbot - ex_excess_n - excess_n, lda);
			utils::copy (1, &data_temp [0] + lda - 2 * nbot - ex_excess_n - excess_n, &data_temp [0] + lda - nbot, lda, lda);
		} else {
			utils::copy (1, &values_n [0], &data_temp [n - 1 + ntop + ex_excess_0], 1, lda);
		}

		utils::matrix_add_scaled (n - 2 + ntop + nbot - excess_0 - excess_n, 1, 1.0, data + 1 - ntop + excess_0, &data_temp [ex_excess_0 + 1 + excess_0], n, lda);
		utils::interpolate (ex_excess_0, 1, n, 1.0, 1.0, positions, data, &positions_0 [0], &data_temp [1], n, lda);
		utils::interpolate (ex_excess_n, 1, n, 1.0, 1.0, positions, data, &positions_n [0], &data_temp [lda - 1 - ex_excess_n], n, lda);
		
		utils::p_block_matrix_solve (messenger_ptr->get_id (), messenger_ptr->get_np (), n - excess_0 - excess_n - ntop - nbot, excess_0 + ex_excess_0 + 2 * ntop, excess_n + ex_excess_n + 2 * nbot, &factorized_matrix [0], &ipiv [0], &data_temp [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, 1, lda, sqrt ((int) boundary_matrix.size ()), lda);
		
		TRACE ("Matrix solve complete.");
		
		for (int i = 0; i < 1; ++i) {
			for (int j = 0; j < n; ++j) {
				if (std::isnan (data_temp [i * n + j])) {
					throw exceptions::nan ();
				}
			}
		}
				
		TRACE ("Updating...");
		utils::matrix_copy (n, 1, &data_temp [ex_excess_0 + ntop], data, lda, n);
		
		*component_flags |= transformed_vertical;
	}

	template class solver <double>;
} /* one_d */