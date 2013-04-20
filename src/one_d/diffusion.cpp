/*!***********************************************************************
 * \file one_d/diffusion.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <iostream>
#include <cmath>
#include <cassert>
#include <memory>
#include "../config.hpp"
#include "diffusion.hpp"
#include "../collocation/collocation.hpp"

namespace one_d
{
	explicit_diffusion::explicit_diffusion (double i_coeff, double *i_timestep_ptr, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_data_in, double *i_data_out, int *i_flags_ptr) : bases::explicit_plan (i_n, i_data_in, i_data_out, i_flags_ptr) 
	{
		coeff = i_coeff;
		timestep_ptr = i_timestep_ptr;
		grid = i_grid;
	}

	void explicit_diffusion::execute () {
		int ione = 1;
		char charN = 'N';
		double dpone = 1.0;
	
		TRACE ("Operating...");
		
		double scalar = coeff * *timestep_ptr;
		// Set up and evaluate the explicit part of the diffusion equation
		dgemv_ (&charN, &n, &n, &scalar, grid->get_data (2), &n, data_in, &ione, &dpone, data_out, &ione);

		TRACE ("Operation complete.");
	}

	implicit_diffusion::implicit_diffusion (double i_coeff, double i_alpha_0, double i_alpha_n, double *i_timestep_ptr, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_matrix, int *i_flags_ptr) : bases::implicit_plan (i_n, i_matrix, i_flags_ptr) {
		coeff = i_coeff;
		alpha_0 = i_alpha_0;
		alpha_n = i_alpha_n;
		timestep_ptr = i_timestep_ptr;
		n = i_n;
		grid = i_grid;
		matrix = i_matrix;
	}

	void implicit_diffusion::execute () {			
		TRACE ("Operating...");

		double scalar = coeff * *timestep_ptr * alpha_0;
	   	daxpy_ (&n, &scalar, grid->get_data (2), &n, matrix, &n);

		scalar = coeff * *timestep_ptr;
		// This is the main loop for setting up the diffusion equation in Chebyshev space
		for (int i = 1; i < n - 1; ++i) {
		   	daxpy_ (&n, &scalar, grid->get_data (2) + i, &n, matrix + i, &n);
		}
		
		scalar = coeff * *timestep_ptr * alpha_n;
	   	daxpy_ (&n, &scalar, grid->get_data (2) + n - 1, &n, matrix + n - 1, &n);
	
		TRACE ("Operation complete.");
	}
} /* one_d */

