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
#include "../bases/element.hpp"
#include "element_one_d.hpp"

namespace one_d
{
	solver::solver (bases::element* i_element_ptr, int i_n, double& i_timestep, double& i_alpha_0, double& i_alpha_n, double *i_default_matrix, double *i_matrix, int i_name_in, int i_name_rhs, int i_name_out) : bases::solver (i_element_ptr, i_n, i_name_in, i_name_out), timestep (i_timestep), alpha_0 (i_alpha_0), alpha_n (i_alpha_n) {
		error.resize (2, 0.0);
		data_temp.resize (n, 0.0);
		rhs = &((*element_ptr) [i_name_rhs]);
		default_matrix = i_default_matrix;
		matrix = i_matrix;
		factorized_matrix.resize (n * n);
		ipiv.resize (n, 0);
	}
	
	void solver::factorize () {
		int info;
		
		bases::solver::factorize ();
		
		utils::copy (n * n, default_matrix, &factorized_matrix [0]);
		
		utils::add_scaled (n, alpha_0 * timestep, matrix + 0, &factorized_matrix [0] + 0, n, n);	
		for (int i = 1; i < n - 1; ++i) {
			utils::add_scaled (n, timestep, matrix + i, &factorized_matrix [0] + i, n, n);	
		}
		utils::add_scaled (n, alpha_n * timestep, matrix + n - 1, &factorized_matrix [0] + n - 1, n, n);

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
		switch (*status) {
			case 0:
				for (int i = 0; i < n; ++i) {
					utils::copy (n, default_matrix + i * n, global_matrix + my_index * (N + 1) + i * N);
				}
				utils::add_scaled (n, alpha_0 * timestep, matrix, global_matrix + my_index * (N + 1), n, N);
				for (int i = 1; i < n - 1; ++i) {
					utils::add_scaled (n, timestep, matrix + i, global_matrix + my_index * (N + 1) + i, n, N);
				}
				utils::add_scaled (n, alpha_n * timestep, matrix + n - 1, global_matrix + my_index * (N + 1) + n - 1, n, N);
				if (element_ptr->is_linked (edge_0)) {
					element_ptr->send (n, matrix, edge_0, alpha_0 * timestep, n);
				}
				if (element_ptr->is_linked (edge_n)) {
					element_ptr->send (n, matrix + n - 1, edge_n, alpha_n * timestep, n);
				}
				*status = -1;
				break;
			case -1:
				if (element_ptr->is_linked (edge_0)) {
					element_ptr->recv (n, global_matrix + N * boundary_0 + my_index, edge_0, 0.0, N);
				}
				if (element_ptr->is_linked (edge_n)) {
					element_ptr->recv (n, global_matrix + N * boundary_n + my_index + n - 1, edge_n, 0.0, N);
				}
				*status = -2;
				break;
			case -2:
				utils::copy (n, data_in, global_rhs + my_index);
				utils::add_scaled (n, timestep, rhs, global_rhs + my_index);
				element_ptr->send (1, global_rhs + my_index, edge_0);
				element_ptr->send (1, global_rhs + my_index + n - 1, edge_n);
				*status = -3;
				break;
			case -3:
				element_ptr->recv (1, global_rhs + my_index, edge_0);
				element_ptr->recv (1, global_rhs + my_index + n - 1, edge_n);
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
		
		MDEBUG ("error " << error [edge_0] << " " << error [edge_n]);
		
		utils::copy (n, data_in, &data_temp [0]);
		data_temp [0] += alpha_0 * timestep * rhs [0] + error [edge_0];
		data_temp [n - 1] += alpha_n * timestep * rhs [n - 1] + error [edge_n];
		utils::add_scaled (n - 2, timestep, rhs + 1, &data_temp [1]);
		
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
	
	void solver::calculate_bounds () {
		data_0 = utils::dot (n, default_matrix, &data_temp [0], n);
		data_n = utils::dot (n, default_matrix + n - 1, &data_temp [0], n);
		prev_data_0 = data_0;
		prev_data_n = data_n;
	}
	
	void solver::send_bounds () {
		element_ptr->send (1, &data_0, edge_0);
		element_ptr->send (1, &data_n, edge_n);
	}
	
	void solver::recv_bounds () {
		element_ptr->recv (1, &data_0, edge_0);
		element_ptr->recv (1, &data_n, edge_n);
	}
	
	void solver::calculate_error () {
		MDEBUG ("data " << data_0 << " " << prev_data_0);
		
		for (int i = 0; i < n; i++) {
			data_temp [i] += matrix [i] * (data_0 - prev_data_0) * timestep;
			data_temp [i] += matrix [(n - 1) * n + i] * (data_n - prev_data_n) * timestep;
		}
		
		error [edge_0] = timestep * (rhs [0] - utils::dot (n, matrix, &data_temp [0], n));
		error [edge_n] = timestep * (rhs [n - 1] - utils::dot (n, matrix + n - 1, &data_temp [0], n));
	}
	
	void solver::send_error () {
		element_ptr->send (1, &error [edge_0], edge_0);
		element_ptr->send (1, &error [edge_n], edge_n);
	}
	
	void solver::recv_error () {
		element_ptr->recv (1, &error [edge_0], edge_0, 0.0);
		element_ptr->recv (1, &error [edge_n], edge_n, 0.0);
	}
	
	void solver::update () {
		error [edge_0] = 0.0;
		error [edge_n] = 0.0;
		utils::copy (n, &data_temp [0], data_out);
	}
} /* one_d */