/*!***********************************************************************
 * \file solver_one_d.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-13.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cmath>
#include "../config.hpp"
#include "solver_one_d.hpp"
#include "../utils/utils.hpp"
#include "../utils/solver_utils.hpp"
#include "../utils/interpolate.hpp"
#include "../bases/element.hpp"
#include "element_one_d.hpp"

namespace one_d
{
	solver::solver (bases::element* i_element_ptr, int i_n, double& i_timestep, double& i_alpha_0, double& i_alpha_n, double *i_default_matrix, double *i_matrix, int i_name_in, int i_name_rhs, int i_name_out) : bases::solver (i_element_ptr, i_n, i_name_in, i_name_out), timestep (i_timestep), alpha_0 (i_alpha_0), alpha_n (i_alpha_n) {
		data_temp.resize (n, 0.0);
		rhs = &((*element_ptr) [i_name_rhs]);
		default_matrix = i_default_matrix;
		matrix = i_matrix;
		factorized_matrix.resize (n * n);
		ipiv.resize (n, 0);
		expected_excess_0 = element_ptr->get_expected_excess (edge_0);
		expected_excess_n = element_ptr->get_expected_excess (edge_n);
		excess_0 = element_ptr->get_excess (edge_0);
		excess_n = element_ptr->get_excess (edge_n);
		if (expected_excess_0 > excess_0) {
			error_0.resize (expected_excess_0 + 1, 0.0);
		} else {
			error_0.resize (excess_0 + 1, 0.0);
		}
		if (expected_excess_n > excess_n) {
			error_n.resize (expected_excess_n + 1, 0.0);
		} else {
			error_n.resize (excess_n + 1, 0.0);
		}
	}
	
	void solver::factorize () {
		int info;
		
		bases::solver::factorize ();
		
		utils::copy (n * n, default_matrix, &factorized_matrix [0]);
		
		utils::add_scaled (n, alpha_0 * timestep, matrix + excess_0, &factorized_matrix [excess_0], n, n);	
		for (int i = excess_0 + 1; i < n - excess_n - 1; ++i) {
			utils::add_scaled (n, timestep, matrix + i, &factorized_matrix [i], n, n);	
		}
		utils::add_scaled (n, alpha_n * timestep, matrix + n - 1 - excess_n, &factorized_matrix [n - 1 - excess_n], n, n);

		utils::matrix_factorize (n, n, &factorized_matrix [0], &ipiv [0], &info);
		
		if (info != 0) {
			ERROR (logger, "Unable to invert matrix");
			throw 0; // For now, kill the program. 
			/*
				TODO Replace this with a more useful exception that can be handled
			*/
		}
	}
	
	void solver::update_globals (int N, double* global_matrix, double* global_rhs, int* status) {
		int my_index = element_ptr->get_index ();
		int boundary_0 = element_ptr->get_boundary_index (edge_0);
		int boundary_n = element_ptr->get_boundary_index (edge_n);
		int excess_0 = element_ptr->get_excess (edge_0);
		int excess_n = element_ptr->get_excess (edge_n);
		std::vector <double> rhs_0 (expected_excess_0, 0.0);
		std::vector <double> rhs_n (expected_excess_n, 0.0);
		std::vector <double> matrix_0 (expected_excess_0 * n, 0.0);
		std::vector <double> matrix_n (expected_excess_n * n, 0.0);
		switch (*status) {
			case 0:
				for (int i = 0; i < n; ++i) {
					if (i < excess_0 || i >= n - excess_n) {
						utils::add_scaled (n, -1.0, default_matrix + i, global_matrix + my_index * (N + 1) + i, n, N);
					} else {
						utils::copy (n, default_matrix + i, global_matrix + my_index * (N + 1) + i, n, N);
					}
				}
				utils::add_scaled (n, alpha_0 * timestep, matrix + excess_0, global_matrix + my_index * (N + 1) + excess_0, n, N);
				for (int i = 1 + excess_0; i < n - 1 - excess_n; ++i) {
					utils::add_scaled (n, timestep, matrix + i, global_matrix + my_index * (N + 1) + i, n, N);
				}
				utils::add_scaled (n, alpha_n * timestep, matrix + n - 1 - excess_n, global_matrix + my_index * (N + 1) + n - 1 - excess_n, n, N);
				if (element_ptr->is_linked (edge_0)) {
					element_ptr->send (n, alpha_0 * timestep, matrix + excess_0, edge_0, n);
				}
				if (element_ptr->is_linked (edge_n)) {
					element_ptr->send (n, alpha_n * timestep, matrix + n - 1 - excess_n, edge_n, n);
				}
				*status = -1;
				break;
			case -1:
				if (element_ptr->is_linked (edge_0)) {
					element_ptr->recv (n, 0.0, global_matrix + N * boundary_0 + my_index + excess_0, edge_0, N);
				}
				if (element_ptr->is_linked (edge_n)) {
					element_ptr->recv (n, 0.0, global_matrix + N * boundary_n + my_index + n - 1 - excess_n, edge_n, N);
				}
				*status = -6;
				break;
			case -6:
				element_ptr->send (1, 1, &excess_0, edge_0);
				element_ptr->send (1, 1, &excess_n, edge_n);
				*status = -7;
				break;
			case -7:
				element_ptr->recv (1, 0, &expected_excess_0, edge_0);
				element_ptr->recv (1, 0, &expected_excess_n, edge_n);
				*status = -8;
				break;
			case -8:
				if (element_ptr->is_linked (edge_0)) {
					element_ptr->send (excess_0, 1.0, &((*element_ptr) (position)), edge_0);
				}
				if (element_ptr->is_linked (edge_n)) {
					element_ptr->send (excess_n, 1.0, &((*element_ptr) (position, n - excess_n)), edge_n);
				}
				*status = -9;
				break;
			case -9:
				if (element_ptr->is_linked (edge_0)) {
					positions_0.resize (expected_excess_0);
					element_ptr->recv (expected_excess_0, 0.0, &positions_0 [0], edge_0);
				}
				if (element_ptr->is_linked (edge_n)) {
					positions_n.resize (expected_excess_n);
					element_ptr->recv (expected_excess_n, 0.0, &positions_n [0], edge_n);
				}
				*status = -12;
				break;
			// case -10:
			// 	for (int i = 0; i < expected_excess_0; ++i) {
			// 		rhs_0 [i] = utils::interpolate (n, &((*element_ptr) (position)), data_in, positions_0 [i]);
			// 		rhs_0 [i] += utils::interpolate (n, &((*element_ptr) (position)), rhs, positions_0 [i]);
			// 	}
			// 	for (int i = 0; i < expected_excess_n; ++i) {
			// 		rhs_n [i] = utils::interpolate (n, &((*element_ptr) (position)), data_in, positions_n [i]);
			// 		rhs_n [i] += utils::interpolate (n, &((*element_ptr) (position)), rhs, positions_n [i]);
			// 	}
			// 	if (element_ptr->is_linked (edge_0)) {
			// 		element_ptr->send (expected_excess_0, 1.0, &rhs_0 [0], edge_0);
			// 	}
			// 	if (element_ptr->is_linked (edge_n)) {
			// 		element_ptr->send (expected_excess_n, 1.0, &rhs_n [0], edge_n);
			// 	}
			// 	*status = -11;
			// 	break;
			// case -11:
			// 	if (element_ptr->is_linked (edge_0)) {
			// 		element_ptr->recv (excess_0, 0.0, global_rhs + my_index, edge_0);
			// 	}
			// 	if (element_ptr->is_linked (edge_n)) {
			// 		element_ptr->recv (excess_n, 0.0, global_rhs + my_index + n - excess_n, edge_n);
			// 	}
			// 	*status = -12;
			// 	break;
			case -12:
				for (int i = 0; i < expected_excess_0; ++i) {
					utils::matrix_interpolate (n, &((*element_ptr) (position)), n, default_matrix, 0.0, expected_excess_0, &matrix_0 [0], positions_0 [i]);
					// utils::matrix_interpolate (n, &((*element_ptr) (position)), n, matrix, 1.0, expected_excess_0, &matrix_0 [0], positions_0 [i]);
				}
				for (int i = 0; i < expected_excess_n; ++i) {
					utils::matrix_interpolate (n, &((*element_ptr) (position)), n, default_matrix, 0.0, expected_excess_n, &matrix_n [0], positions_n [i]);
					// utils::matrix_interpolate (n, &((*element_ptr) (position)), n, matrix, 1.0, expected_excess_n, &matrix_n [0], positions_n [i]);
				}
				element_ptr->send (expected_excess_0 * n, 1.0, &matrix_0 [0], edge_0);
				element_ptr->send (expected_excess_n * n, 1.0, &matrix_n [0], edge_n);
				*status = -13;
				break;
			case -13:
				if (element_ptr->is_linked (edge_0)) {
					element_ptr->recv (excess_0 * n, 0.0, &matrix_0 [0], edge_0);
				}
				if (element_ptr->is_linked (edge_n)) {
					element_ptr->recv (excess_n * n, 0.0, &matrix_n [0], edge_n);
				}
				for (int i = 0; i < excess_0; ++i) {
					utils::copy (n, &matrix_0 [i * n], global_matrix + my_index + i + boundary_0 * N, 1, N);
				}
				for (int i = 0; i < excess_n; ++i) {
					utils::copy (n, &matrix_n [i * n], global_matrix + my_index + n - excess_n + i + boundary_n * N, 1, N);
				}
				*status = -2;
				/*
					TODO Adjacent elements must have the same number of indices, make more general
				*/
			case -2:
				utils::copy (n - excess_0 - excess_n, data_in + excess_0, global_rhs + my_index + excess_0);
				utils::add_scaled (n - excess_0 - excess_n, timestep, rhs + excess_0, global_rhs + my_index + excess_0);
				element_ptr->send (1, global_rhs + my_index + excess_0, edge_0);
				element_ptr->send (1, global_rhs + my_index + n - 1 - excess_n, edge_n);
				*status = -3;
				break;
			case -3:
				element_ptr->recv (1, global_rhs + my_index + excess_0, edge_0);
				element_ptr->recv (1, global_rhs + my_index + n - 1 - excess_n, edge_n);
				*status = 1;
				break;
			default:
				break;
		}
	}
	
	void solver::update_from_globals (double* global_out) {
		utils::copy (n, global_out + element_ptr->get_index (), data_out);
	}
		
	void solver::execute () {
		int info;
		
		bases::solver::execute ();
		
		TRACE (logger, "Solving...");
		
		MDEBUG ("error " << error_0 [0] << " " << error_n [0]);
		
		utils::copy (n, data_in, &data_temp [0]);
		
		data_temp [excess_0] += alpha_0 * timestep * rhs [excess_0] + error_0 [0];
		for (int i = 0; i < excess_0; ++i) {
			data_temp [i] += error_0 [excess_0 - i]; 
		}
		data_temp [n - 1 - excess_n] += alpha_n * timestep * rhs [n - 1 - excess_n] + error_n [0];
		for (int i = 0; i < excess_n; ++i) {
			data_temp [n - excess_n + i] += error_n [i + 1];
		}
		utils::add_scaled (n - 2 - excess_0 - excess_n, timestep, rhs + 1 + excess_0, &data_temp [excess_0 + 1]);
		
		utils::matrix_solve (n, &factorized_matrix [0], &ipiv [0], &data_temp [0], &info);
		
		for (int i = 0; i < n; ++i) {
			if (std::isnan (data_temp [i])) {
				info = -1;
			}
		}
		
		if (info != 0) {
			ERROR (logger, "Unable to solve factorized matrix equation");
			throw 0; // For now, kill the program. 
			/*
				TODO Replace this with a more useful exception that can be handled
			*/
		}
		TRACE (logger, "Solve complete.")
	}
	
	void solver::send_positions () {
		TRACE (logger, "Sending positions...");
		
		if (excess_0 != 0 && (int) positions_0.size () == 0) {
			element_ptr->send (excess_0, 1.0, &((*element_ptr) (position)), edge_0);
		} 
		if (excess_n != 0 && (int) positions_n.size () == 0) {
			element_ptr->send (excess_n, 1.0, &((*element_ptr) (position, n - 1)), edge_n, -1);
		}
	}
	
	void solver::recv_positions () {
		TRACE (logger, "Recving positions...");
		
		if (expected_excess_0 != 0 && (int) positions_0.size () == 0) {
			positions_0.resize (expected_excess_0);
			element_ptr->recv (expected_excess_0, &positions_0 [0], edge_0);
		}
		if (expected_excess_n != 0 && (int) positions_n.size () == 0) {
			positions_n.resize (expected_excess_n);
			element_ptr->recv (expected_excess_n, &positions_n [0], edge_n);
		}
	}
	
	void solver::calculate_bounds () {
		TRACE (logger, "Calculating bounds...");
		
		error_0 [0] = (alpha_0 * timestep) * (rhs [excess_0] - utils::dot (n, matrix + excess_0, &data_temp [0], n));
		error_n [0] = (alpha_n * timestep) * (rhs [n - 1 - excess_n] - utils::dot (n, matrix + n - 1 - excess_n, &data_temp [0], n));
		for (int i = 0; i < excess_0; ++i) {
			error_0 [i + 1] = utils::dot_interpolate (n, &((*element_ptr) (position)), n, default_matrix, &data_temp [0], positions_0 [i]);
			MDEBUG ("CALC_0 " << error_0 [i + 1]);
		}
		for (int i = 0; i < excess_n; ++i) {
			error_n [i + 1] = utils::dot_interpolate (n, &((*element_ptr) (position)), n, default_matrix, &data_temp [0], positions_n [i]);
			MDEBUG ("CALC_N " << error_n [i + 1]);
		}
	}
	
	void solver::send_bounds () {
		TRACE (logger, "Sending bounds...");
		
		element_ptr->send (expected_excess_0 + 1, 1.0, &error_0 [0], edge_0);
		element_ptr->send (expected_excess_n + 1, 1.0, &error_n [0], edge_n);
	}
	
	void solver::recv_bounds () {
		TRACE (logger, "Recving bounds...");
		
		element_ptr->recv (excess_0 + 1, 0.0, &error_0 [0], edge_0);
		element_ptr->recv (excess_n + 1, 0.0, &error_n [0], edge_n);
		MDEBUG ("ERROR_0 " << error_0 [0]);
		for (int i = 0; i < excess_0; ++i) {
			MDEBUG ("before " << error_0 [i + 1]);
			error_0 [i + 1] -= data_in [excess_0 - 1 - i];
			MDEBUG (error_0 [i + 1]);
		}
		MDEBUG ("ERROR_N " << error_n [0]);
		for (int i = 0; i < excess_n; ++i) {
			error_n [i + 1] -= data_in [n - excess_n + i];
			MDEBUG (error_n [i + 1]);
		}
	}
	
	void solver::calculate_error () {
		// MDEBUG ("data " << data_0 << " " << prev_data_0);
		// 
		// for (int i = 0; i < n; i++) {
		// 	data_temp [i] += matrix [i] * (data_0 - prev_data_0) * timestep;
		// 	data_temp [i] += matrix [(n - 1) * n + i] * (data_n - prev_data_n) * timestep;
		// }
		// 
		// error [edge_0] = timestep * (rhs [0] - utils::dot (n, matrix, &data_temp [0], n));
		// error [edge_n] = timestep * (rhs [n - 1] - utils::dot (n, matrix + n - 1, &data_temp [0], n));
	}
	
	void solver::send_error () {
		// element_ptr->send (1, &error [edge_0], edge_0);
		// element_ptr->send (1, &error [edge_n], edge_n);
	}
	
	void solver::recv_error () {
		// element_ptr->recv (1, &error [edge_0], edge_0, 0.0);
		// element_ptr->recv (1, &error [edge_n], edge_n, 0.0);
	}
	
	void solver::update () {
		TRACE (logger, "Updating...");
		// error [edge_0] = 0.0;
		// error [edge_n] = 0.0;
		utils::copy (n, &data_temp [0], data_out);
	}
} /* one_d */