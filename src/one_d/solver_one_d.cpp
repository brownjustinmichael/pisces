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
	template <class datatype>
	solver <datatype>::solver (bases::messenger* i_messenger_ptr, int i_n, int i_excess_0, int i_excess_n, int i_n_iterations, datatype& i_timestep, datatype& i_alpha_0, datatype& i_alpha_n, datatype* i_positions, datatype *i_default_matrix, datatype *i_matrix, datatype* i_data_in, datatype* i_explicit_rhs, datatype* i_implicit_rhs, datatype* i_data_out, int i_flags) : 
	bases::solver <datatype> (i_flags), 
	explicit_plan <datatype> (i_n, i_data_in, i_data_out),
	messenger_ptr (i_messenger_ptr),
	timestep (i_timestep), 
	alpha_0 (i_alpha_0), 
	alpha_n (i_alpha_n), 
	positions (i_positions),
	n_iterations (i_n_iterations),
	excess_0 (i_excess_0), 
	excess_n (i_excess_n) { 
		
		TRACE ("Instantiating...");
		
		expected_excess_0 = 0;
		expected_excess_n = 0;
		explicit_rhs = i_explicit_rhs;
		implicit_rhs = i_implicit_rhs;
		default_matrix = i_default_matrix;
		matrix = i_matrix;		
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
	void solver <datatype>::_factorize () {
		int info;
		
		DEBUG ("factoring..." << matrix [126]);
		
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
	void solver <datatype>::execute () {
		int info;
		int beta_0 = 1.5, beta_1 = -0.5;
		
		bases::solver <datatype>::execute ();
		
		DEBUG ("solving..." << matrix [126]);
		
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
		
		TRACE ("Solve complete.")
	}
	
	template class solver <double>;
	template class solver <float>;
} /* one_d */