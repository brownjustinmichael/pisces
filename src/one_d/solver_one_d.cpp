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

namespace one_d
{
	solver::solver (bases::element* i_element_ptr, int i_n, double& i_timestep, double *i_default_matrix, double *i_matrix, int i_name_in, int i_name_rhs, int i_name_out) : bases::solver (i_element_ptr, i_n, i_name_in, i_name_out), timestep (i_timestep) {
		rhs = &((*element_ptr) [i_name_rhs]);
		default_matrix = i_default_matrix;
		matrix = i_matrix;
		ipiv.resize (n, 0);
	}
	
	void solver::factorize () {
		int info;
		
		bases::solver::factorize ();
		
		utils::scale (n * n, timestep, matrix);
		
		utils::add_scaled (n * n, 1.0, default_matrix, matrix);
		
		utils::matrix_factorize (n, n, matrix, &ipiv [0], &info);
		
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
		
		TRACE (logger, "Solving...")
		
		if (data_out == rhs) {
			utils::add_scaled (n, timestep, data_in, data_out);
		} else {
			if (data_in != data_out) {
				utils::copy (n, data_in, data_out);
			}
			utils::add_scaled (n, timestep, rhs, data_out);
		}
		
		utils::matrix_solve (n, &matrix [0], &ipiv [0], data_out, &info);
		
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
} /* one_d */