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
	solver::solver (bases::element* i_element_ptr, int i_n, int i_excess_0, int i_excess_n, double& i_timestep, double& i_alpha_0, double& i_alpha_n, double *i_default_matrix, double *i_matrix, int i_name_in, int i_name_rhs, int i_name_out) : 
	bases::solver (i_element_ptr, i_n, i_name_in, i_name_out), 
	timestep (i_timestep), 
	alpha_0 (i_alpha_0), 
	alpha_n (i_alpha_n), 
	excess_0 (i_excess_0), 
	excess_n (i_excess_n) { 
		
		TRACE ("Instantiating...");
		
		expected_excess_0 = 0;
		expected_excess_n = 0;
		rhs = &((*element_ptr) [i_name_rhs]);
		default_matrix = i_default_matrix;
		matrix = i_matrix;		
		data_temp.resize (n, 0.0);
		factorized_matrix.resize (n * n);
		ipiv.resize (n, 0);

		TRACE ("Instantiated");
		
		/*
			TODO It would be nice if the boundaries were established by the time this was called so that positions and excesses could be exchanged in the constructor
		*/
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
			ERROR ("Unable to invert matrix");
			throw 0; // For now, kill the program. 
			/*
				TODO Replace this with a more useful exception that can be handled
			*/
		}
	}
	
	void solver::execute () {
		int info;
		
		bases::solver::execute ();
		
		TRACE ("Solving...");
		
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
			ERROR ("Unable to solve factorized matrix equation");
			throw 0; // For now, kill the program. 
			/*
				TODO Replace this with a more useful exception that can be handled
			*/
		}
				
		out_error_0 [0] = (alpha_0 * timestep) * (rhs [excess_0] - utils::dot (n, matrix + excess_0, &data_temp [0], n));
		out_error_n [0] = (alpha_n * timestep) * (rhs [n - 1 - excess_n] - utils::dot (n, matrix + n - 1 - excess_n, &data_temp [0], n));
		for (int i = 0; i < expected_excess_0; ++i) {
			out_error_0 [i + 1] = utils::dot_interpolate (n, &((*element_ptr) (position)), n, default_matrix, &data_temp [0], positions_0 [i]);
		}
		for (int i = 0; i < expected_excess_n; ++i) {
			out_error_n [i + 1] = utils::dot_interpolate (n, &((*element_ptr) (position)), n, default_matrix, &data_temp [0], positions_n [i]);
		}
				
		element_ptr->communicate <double> (expected_excess_0 + 1, edge_0, &out_error_0 [0], expected_excess_n + 1, edge_n, &out_error_n [0], &error_0 [0], &error_n [0], excess_0 + 1, excess_n + 1);
		
		for (int i = 0; i < excess_0; ++i) {
			error_0 [i + 1] -= data_in [i];
		}
		for (int i = 0; i < excess_n; ++i) {
			error_n [i + 1] -= data_in [n - excess_n + i];
		}
		
		TRACE ("Solve complete.")
	}
	
	void solver::send_positions () {
		TRACE ("Sending positions...");
		
		element_ptr->communicate <int> (1, edge_0, &excess_0, 1, edge_n, &excess_n, &expected_excess_0, &expected_excess_n);
		
		error_0.resize (excess_0 + 1, 0.0);
		error_n.resize (excess_n + 1, 0.0);
		out_error_0.resize (expected_excess_0 + 1, 0.0);
		out_error_n.resize (expected_excess_n + 1, 0.0);
		positions_0.resize (expected_excess_0);
		positions_n.resize (expected_excess_n);
		
		element_ptr->communicate <double> (excess_0, edge_0, &((*element_ptr) (position)), excess_n, edge_n, &((*element_ptr) (position, n - 1)), &positions_0 [0], &positions_n [0], expected_excess_0, expected_excess_n);
	}
	
	void solver::update () {
		TRACE ("Updating...");
		utils::copy (n, &data_temp [0], data_out);
	}
} /* one_d */