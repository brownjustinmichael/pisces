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

	void solver::execute () {
		int info;
		
		bases::solver::execute ();
		
		TRACE (logger, "Solving...");
		
		if (data_out == rhs) {
			utils::add_scaled (n, timestep, data_in, data_out);
		} else {
			if (data_in != data_out) {
				utils::copy (n, data_in, data_out);
			}
			data_out [0] += alpha_0 * timestep * rhs [0] + error [edge_0];
			data_out [n - 1] += alpha_n * timestep * rhs [n - 1] + error [edge_n];
			utils::add_scaled (n - 2, timestep, rhs + 1, data_out + 1);
		}
		
		utils::matrix_solve (n, &factorized_matrix [0], &ipiv [0], data_out, &info);
		
		for (int i = 0; i < n; ++i) {
			if (std::isnan (data_out [i])) {
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
	
	void solver::calculate_error () {
/*
		double data_0, data_n;
		double prev_data_0, prev_data_n;
		
		data_0 = utils::dot (n, default_matrix, data_out, n);
		data_n = utils::dot (n, default_matrix + n - 1, data_out, n);
		prev_data_0 = data_0;
		prev_data_n = data_n;
		
		MDEBUG ("Before: " << data_0 << " " << data_n);
		
		element_ptr->send (1, &data_0, edge_0);
		element_ptr->send (1, &data_n, edge_n);
		
		element_ptr->recv (1, &data_0, edge_0);
		element_ptr->recv (1, &data_n, edge_n);
		
		MDEBUG ("After: " << data_0 << " " << data_n);
		
		std::vector <double> temp (n, 0);
		
		data_out [i] += matrix [i] * (boundary_value - boundary_value_0);
		
		new_bit = utils::dot (n, matrix, data_out, n) - rhs [0] * alpha_0;*/

		
		// Send new_bit to adjacent element and subtract it from rhs
		
		// Recalculate data_out using new rhs
		
		// Possible: Add new_bit top and bottom arguments to solve so sending can be handled externally
	}
	
	void solver::update () {
		utils::copy (n, data_out, data_in);
	}
} /* one_d */