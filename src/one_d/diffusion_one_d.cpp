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
#include "diffusion_one_d.hpp"
#include "../utilities/chebyshev.hpp"

namespace one_d
{
	namespace chebyshev
	{
		void explicit_diffusion::init (double i_coeff, double i_alpha_0, double i_alpha_n, double& i_timestep, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_data_in, double *i_data_out, int *i_flags_ptr) {
			coeff = i_coeff;
			alpha_0 = i_alpha_0;
			alpha_n = i_alpha_n;
			grid = i_grid;
		}

		void explicit_diffusion::execute () {
			int ione = 1, nm2 = n - 2;
			char charN = 'N';
			double dpone = 1.0;
		
			TRACE (logger, "Operating...");
		
			double scalar = coeff * timestep * alpha_0;
			data_out [0] += scalar * ddot_ (&n, grid->get_data (2), &n, data_in, &ione);
		
			scalar = coeff * timestep;
			// Set up and evaluate the explicit part of the diffusion equation
			dgemv_ (&charN, &nm2, &n, &scalar, grid->get_data (2) + 1, &n, data_in, &ione, &dpone, data_out + 1, &ione);

			scalar = coeff * timestep * alpha_n;
			data_out [n - 1] += scalar * ddot_ (&n, grid->get_data (2) + n - 1, &n, data_in, &ione);

			TRACE (logger, "Operation complete.");
		}

		implicit_diffusion::implicit_diffusion (double i_coeff, double i_alpha_0, double i_alpha_n, double& i_timestep, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_matrix, int *i_flags_ptr, int i_logger) : bases::implicit_plan (i_n, i_matrix, i_flags_ptr, i_logger), timestep (i_timestep) {
			coeff = i_coeff;
			alpha_0 = i_alpha_0;
			alpha_n = i_alpha_n;
			timestep = i_timestep;
			n = i_n;
			grid = i_grid;
			matrix = i_matrix;
		}

		void implicit_diffusion::execute () {			
			TRACE (logger, "Operating...");

			double scalar = coeff * timestep * alpha_0;
		
		   	daxpy_ (&n, &scalar, grid->get_data (2), &n, matrix, &n);

			scalar = coeff * timestep;
			// This is the main loop for setting up the diffusion equation in Chebyshev space
			for (int i = 1; i < n - 1; ++i) {
			   	daxpy_ (&n, &scalar, grid->get_data (2) + i, &n, matrix + i, &n);
			}
		
			scalar = coeff * timestep * alpha_n;
		   	daxpy_ (&n, &scalar, grid->get_data (2) + n - 1, &n, matrix + n - 1, &n);
	
			TRACE (logger, "Operation complete.");
		}
	} /* chebyshev */
} /* one_d */

