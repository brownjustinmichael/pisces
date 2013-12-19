/*!***********************************************************************
 * \file solver_one_d.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-13.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "solver_one_d.hpp"

#include "../bases/messenger.hpp"
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
	solver <datatype>::solver (bases::grid <datatype> &i_grid, bases::messenger* i_messenger_ptr, datatype& i_timestep, datatype& i_alpha_0, datatype& i_alpha_n, datatype* i_data_in, datatype* i_explicit_rhs, datatype* i_implicit_rhs, datatype* i_data_out, int i_flags) :
	bases::solver <datatype> (i_flags), 
	explicit_plan <datatype> (i_grid, i_data_in, i_data_out),
	messenger_ptr (i_messenger_ptr),
	timestep (i_timestep), 
	alpha_0 (i_alpha_0), 
	alpha_n (i_alpha_n), 
	positions (&(grid.position ())),
	excess_0 (grid.excess_0), 
	excess_n (grid.excess_n),
	explicit_rhs (i_explicit_rhs),
	implicit_rhs (i_implicit_rhs),
	default_matrix (i_grid.get_data (0)) {
		TRACE ("Instantiating...");
		matrix.resize (n * n);
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
		ipiv.resize (n); // Should be n - ntop - nbot - excess_0 - excess_n
		data_temp.resize (n + ex_excess_0 + ex_excess_n + nbot + ntop);
		TRACE ("Instantiated.");
	}
	
	template <class datatype>
	void solver <datatype>::_factorize () {
		int info, lda = n + ex_excess_0 + ex_excess_n + nbot + ntop;
		TRACE ("Factorizing..." << messenger_ptr->get_id ());
		std::stringstream debug;
		
		utils::matrix_scale (lda, lda, 0.0, &factorized_matrix [0], lda);
		utils::matrix_copy (n, n, default_matrix, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1)], n, lda);

		utils::matrix_add_scaled (n - excess_n - excess_0 - 2, n, timestep, &matrix [excess_0 + 1], &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + 1 + excess_0], n, lda);
		if (ntop != 0) {
			utils::matrix_add_scaled (ntop, n, alpha_0 * timestep, &matrix [excess_0], &factorized_matrix [(ntop + ex_excess_0) * lda], n, n + ex_excess_0 + ex_excess_n + ntop + nbot);
			utils::interpolate (ex_excess_0, n, n, timestep, positions, &matrix [0], &positions_0 [0], &factorized_matrix [(ntop + ex_excess_0) * lda + ntop], n, lda);
			utils::matrix_add_scaled (ntop, n, alpha_0 * timestep, &matrix [excess_0], &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + excess_0], n, n + ex_excess_0 + ex_excess_n + ntop + nbot);
		}
		if (nbot != 0) {
			utils::matrix_add_scaled (nbot, n, alpha_n * timestep, &matrix [n - nbot - excess_n], &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + n - nbot - excess_n], n, lda);
			utils::interpolate (ex_excess_n, n, n, timestep, positions, &matrix [0], &positions_n [0], &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + n], n, lda);
			utils::matrix_add_scaled (nbot, n, alpha_n * timestep, &matrix [n - nbot - excess_n], &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + n + ex_excess_n], n, lda);
		}
		

		
		utils::p_block_matrix_factorize (messenger_ptr->get_id (), messenger_ptr->get_np (), n - excess_0 - excess_n - ntop - nbot, excess_0 + ex_excess_0 + 2 * ntop, excess_n + ex_excess_n + 2 * nbot, &factorized_matrix [0], &ipiv [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, lda, sqrt ((int) boundary_matrix.size ()));
		
		if (info != 0) {
			ERROR ("Unable to invert matrix");
			throw 0; // For now, kill the program. 
			/*
				TODO Replace this with a more useful exception that can be handled
			*/
		}
		TRACE ("Done.");
	}
	
	template <class datatype>
	void solver <datatype>::execute (int &element_flags) {
		int info, lda = n + ex_excess_0 + ex_excess_n + nbot + ntop;
		TRACE ("Executing solve...");
		utils::scale (lda, 0.0, &data_temp [0]);
		
		if (!(flags & first_run)) {
			value_0 = data_in [0];
			value_n = data_in [n - 1];
			flags |= first_run;
		}
		
		utils::add_scaled (n - excess_0 - excess_n, timestep, implicit_rhs + excess_0, &data_temp [ex_excess_0 + ntop + excess_0]);
		utils::add_scaled (n - excess_0 - excess_n, timestep, explicit_rhs + excess_0, &data_temp [ex_excess_0 + ntop + excess_0]);
		utils::interpolate (ex_excess_0, 1, n - excess_0 - excess_n, 1.0, positions + excess_0, &data_temp [ex_excess_0 + ntop + excess_0], &positions_0 [0], &data_temp [ntop]);
		utils::interpolate (ex_excess_n, 1, n - excess_0 - excess_n, 1.0, positions + excess_0, &data_temp [ex_excess_0 + ntop + excess_0], &positions_n [0], &data_temp [lda - nbot - ex_excess_n]);
		if (ntop != 0) {
			data_temp [ntop + ex_excess_0 + excess_0] *= alpha_0;
			data_temp [0] = data_temp [ntop + ex_excess_0 + excess_0];
			// data_temp [ntop + ex_excess_0 + excess_0] += data_in [excess_0];
		} else {
			data_temp [0] = value_0;
		}
		if (nbot != 0) {
			data_temp [lda - 1 - nbot - ex_excess_n - excess_n] *= alpha_n;
			data_temp [lda - 1] = data_temp [lda - 1 - nbot - ex_excess_n - excess_n];
			// data_temp [lda - 1 - nbot - ex_excess_n - excess_n] += data_in [n - 1 - excess_n];
		} else {
			data_temp [n - 1 + ntop + ex_excess_0] = value_n;
		}

		utils::add_scaled (n - 2 + ntop + nbot - excess_0 - excess_n, 1.0, data_in + 1 - ntop + excess_0, &data_temp [ex_excess_0 + 1 + excess_0]);
		utils::interpolate (ex_excess_0, 1, n, 1.0, positions, data_in, &positions_0 [0], &data_temp [1]);
		utils::interpolate (ex_excess_n, 1, n, 1.0, positions, data_in, &positions_n [0], &data_temp [lda - 1 - ex_excess_n]);

		utils::p_block_matrix_solve (messenger_ptr->get_id (), messenger_ptr->get_np (), n - excess_0 - excess_n - ntop - nbot, excess_0 + ex_excess_0 + 2 * ntop, excess_n + ex_excess_n + 2 * nbot, &factorized_matrix [0], &ipiv [0], &data_temp [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, 1, lda, sqrt ((int) boundary_matrix.size ()));
		
		utils::copy (n, &data_temp [ex_excess_0 + ntop], data_out);
		
		element_flags |= transformed_vertical;
	}
	
	template class solver <double>;
} /* one_d */