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
	iterative_solver <datatype>::iterative_solver (bases::grid <datatype> &i_grid, bases::messenger* i_messenger_ptr, int i_n_iterations, datatype& i_timestep, datatype& i_alpha_0, datatype& i_alpha_n, datatype* i_data_in, datatype* i_explicit_rhs, datatype* i_implicit_rhs, datatype* i_data_out, int i_flags) : 
	bases::solver <datatype> (i_flags), 
	explicit_plan <datatype> (i_grid, i_data_in, i_data_out),
	messenger_ptr (i_messenger_ptr),
	timestep (i_timestep), 
	alpha_0 (i_alpha_0), 
	alpha_n (i_alpha_n), 
	positions (&(grid.position ())),
	n_iterations (i_n_iterations),
	excess_0 (grid.excess_0), 
	excess_n (grid.excess_n) { 
		
		TRACE ("Instantiating...");
		
		expected_excess_0 = 0;
		expected_excess_n = 0;
		explicit_rhs = i_explicit_rhs;
		implicit_rhs = i_implicit_rhs;
		default_matrix = grid.get_data (0);
		matrix = grid.matrix_ptr ();		
		data_temp.resize (n, 0.0);
		factorized_matrix.resize (n * n);
		ipiv.resize (n, 0);
		previous_rhs.resize (n, 0);
		
		messenger_ptr->template send <int> (1, &excess_0, edge_0);
		messenger_ptr->template send <int> (1, &excess_n, edge_n);
		messenger_ptr->template recv <int> (1, &expected_excess_0, edge_0);
		messenger_ptr->template recv <int> (1, &expected_excess_n, edge_n);
		
		error_0.resize (excess_0 + 1, 0.0);
		error_n.resize (excess_n + 1, 0.0);
		out_error_0.resize (expected_excess_0 + 1, 0.0);
		out_error_n.resize (expected_excess_n + 1, 0.0);
		positions_0.resize (expected_excess_0);
		positions_n.resize (expected_excess_n);
		
		messenger_ptr->template send <datatype> (excess_0, positions, edge_0);
		messenger_ptr->template send <datatype> (excess_n, &(positions [n - 1]), edge_n);
		messenger_ptr->template recv <datatype> (expected_excess_0, &positions_0 [0], edge_0);
		messenger_ptr->template recv <datatype> (expected_excess_n, &positions_n [0], edge_n);
		
		TRACE ("Instantiated");
	}
	
	template <class datatype>
	void iterative_solver <datatype>::_factorize () {
		int info;
		std::stringstream debug;
		
		TRACE ("Factorizing...");
		
		utils::copy (n * n, default_matrix, &factorized_matrix [0]);

		utils::add_scaled (n, alpha_0 * timestep, matrix + excess_0, &factorized_matrix [excess_0], n, n);	
		for (int i = excess_0 + 1; i < n - excess_n - 1; ++i) {
			utils::add_scaled (n, timestep, matrix + i, &factorized_matrix [i], n, n);	
		}
		utils::add_scaled (n, alpha_n * timestep, matrix + n - 1 - excess_n, &factorized_matrix [n - 1 - excess_n], n, n);

		utils::matrix_factorize (n, n, &factorized_matrix [0], &ipiv [0], &info);
		
		if (info != 0) {
			ERROR ("Unable to invert matrix");
			throw 0; // For now, kill the program. 
			/*
				TODO Replace this with a more useful exception that can be handled
			*/
		}
	}
	
	template <class datatype>
	void iterative_solver <datatype>::execute (int &element_flags) {
		int info;
		double beta_0 = 1.5, beta_1 = -0.5;
		
		bases::solver <datatype>::execute (element_flags);
		
		TRACE ("Solving...");
		
		for (int j = 0; j < n_iterations; ++j) {
			if (j != 0) {
				if (flags & first_run) {
					out_error_0 [0] = (alpha_0 * timestep) * (explicit_rhs [excess_0] + implicit_rhs [excess_0] - utils::dot (n, matrix + excess_0, &data_temp [0], n));
					out_error_n [0] = (alpha_n * timestep) * (explicit_rhs [n - 1 - excess_n] + implicit_rhs [n - 1 - excess_n] - utils::dot (n, matrix + n - 1 - excess_n, &data_temp [0], n));
				} else {
					out_error_0 [0] = (alpha_0 * timestep) * (beta_0 * explicit_rhs [excess_0] + beta_1 * previous_rhs [excess_0] + implicit_rhs [excess_0] - utils::dot (n, matrix + excess_0, &data_temp [0], n));
					out_error_n [0] = (alpha_n * timestep) * (beta_0 * explicit_rhs [n - 1 - excess_n] + beta_1 * previous_rhs [n - 1 - excess_n] + implicit_rhs [n - 1 - excess_n] - utils::dot (n, matrix + n - 1 - excess_n, &data_temp [0], n));
				}
			
				for (int i = 0; i < expected_excess_0; ++i) {
					out_error_0 [i + 1] = utils::dot_interpolate (n, positions, n, default_matrix, &data_temp [0], positions_0 [i]);
				}
				for (int i = 0; i < expected_excess_n; ++i) {
					out_error_n [i + 1] = utils::dot_interpolate (n, positions, n, default_matrix, &data_temp [0], positions_n [i]);
				}
					
				messenger_ptr->template send <datatype> (expected_excess_0 + 1, &out_error_0 [0], edge_0);
				messenger_ptr->template send <datatype> (expected_excess_n + 1, &out_error_n [0], edge_n);
				messenger_ptr->template recv <datatype> (excess_0 + 1, &error_0 [0], edge_0);
				messenger_ptr->template recv <datatype> (excess_n + 1, &error_n [0], edge_n);
					
				for (int i = 0; i < excess_0; ++i) {
					error_0 [i + 1] -= data_in [i];
				}
				for (int i = 0; i < excess_n; ++i) {
					error_n [i + 1] -= data_in [n - excess_n + i];
				}	
			}
			
			utils::copy (n, data_in, &data_temp [0]);
		
			if (flags & first_run) {
				data_temp [excess_0] += alpha_0 * timestep * (explicit_rhs [excess_0] + implicit_rhs [excess_0]) + error_0 [0];
				data_temp [n - 1 - excess_n] += alpha_n * timestep * (explicit_rhs [n - 1 - excess_n] + implicit_rhs [n - 1 - excess_n]) + error_n [0];
				utils::add_scaled (n - 2 - excess_0 - excess_n, timestep, explicit_rhs + 1 + excess_0, &data_temp [excess_0 + 1]);
			} else {
				data_temp [excess_0] += alpha_0 * timestep * (beta_0 * explicit_rhs [excess_0] + beta_1 * previous_rhs [excess_0] + implicit_rhs [excess_0]) + error_0 [0];
				data_temp [n - 1 - excess_n] += alpha_n * timestep * (beta_0 * explicit_rhs [n - 1 - excess_n] + beta_1 * explicit_rhs [n - 1 - excess_n] + beta_0 * implicit_rhs [n - 1 - excess_n]) + error_n [0];
				utils::add_scaled (n - 2 - excess_0 - excess_n, beta_0 * timestep, explicit_rhs + 1 + excess_0, &data_temp [excess_0 + 1]);
				utils::add_scaled (n - 2 - excess_0 - excess_n, beta_1 * timestep, &previous_rhs [1 + excess_0], &data_temp [excess_0 + 1]);
			}
			
			for (int i = 0; i < excess_0; ++i) {
				data_temp [i] += error_0 [excess_0 - i]; 
			}
			for (int i = 0; i < excess_n; ++i) {
				data_temp [n - excess_n + i] += error_n [i + 1];
			}
			utils::add_scaled (n - 2 - excess_0 - excess_n, timestep, implicit_rhs + 1 + excess_0, &data_temp [excess_0 + 1]);

			utils::matrix_solve (n, &factorized_matrix [0], &ipiv [0], &data_temp [0], &info);
		
			for (int i = 0; i < n; ++i) {
				if (std::isnan (data_temp [i])) {
					ERROR ("Nan detected.");
					info = -1;
				}
			}
		
			if (info != 0) {
				ERROR ("Unable to solve factorized matrix equation");
				throw 0; // For now, kill the program. 
				/*
					TODO Replace this with a more useful exception that can be handled
				*/
			}
			
		}
		TRACE ("Updating...");
		utils::copy (n, &data_temp [0], data_out);
		utils::copy (n, explicit_rhs, &previous_rhs [0]);
		flags |= first_run;
		element_flags |= transformed_vertical;
		
		TRACE ("Solve complete.")
	}
	
	template class iterative_solver <double>;
	template class iterative_solver <float>;
	
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
	default_matrix (i_grid.get_data (0)), 
	matrix (i_grid.matrix_ptr ()) {
		TRACE ("Instantiating...");
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
	}
	
	template <class datatype>
	void solver <datatype>::_factorize () {
		int info, lda = n + ex_excess_0 + ex_excess_n + nbot + ntop;
		TRACE ("Factorizing..." << messenger_ptr->get_id ());
		std::stringstream debug;
		
		utils::matrix_copy (n, n, default_matrix, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1)], n, lda);

		utils::matrix_add_scaled (n - excess_n - excess_0 - 2, n, timestep, matrix + excess_0 + 1, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + 1 + excess_0], n, lda);
		if (ntop != 0) {
			utils::matrix_add_scaled (ntop, n, alpha_0 * timestep, matrix + excess_0, &factorized_matrix [(ntop + ex_excess_0) * lda], n, n + ex_excess_0 + ex_excess_n + ntop + nbot);
			utils::interpolate (ex_excess_0, n, n, timestep, positions, matrix, &positions_0 [0], &factorized_matrix [(ntop + ex_excess_0) * lda + ntop], n, lda);
			utils::matrix_add_scaled (ntop, n, alpha_0 * timestep, matrix + excess_0, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + excess_0], n, n + ex_excess_0 + ex_excess_n + ntop + nbot);
		}
		if (nbot != 0) {
			utils::matrix_add_scaled (nbot, n, alpha_n * timestep, matrix + n - nbot - excess_n, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + n - nbot - excess_n], n, lda);
			utils::interpolate (ex_excess_n, n, n, timestep, positions, matrix, &positions_n [0], &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + n], n, lda);
			utils::matrix_add_scaled (nbot, n, alpha_n * timestep, matrix + n - nbot - excess_n, &factorized_matrix [(ntop + ex_excess_0) * (lda + 1) + n + ex_excess_n], n, lda);
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