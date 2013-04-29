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
		void explicit_diffusion::init (double i_coeff, double& i_timestep, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_data_in, double *i_data_out, int *i_flags_ptr) {
			coeff = i_coeff;
			grid = i_grid;
		}

		void explicit_diffusion::execute () {
			int ione = 1, nm2 = n - 2;
			char charN = 'N';
			double dpone = 1.0;
			double scalar = coeff * timestep;
			
			bases::explicit_plan::execute ();
		
			TRACE (logger, "Operating...");
		
			data_out [0] += scalar * ddot_ (&n, grid->get_data (2), &n, data_in, &ione);
		
			// Set up and evaluate the explicit part of the diffusion equation
			dgemv_ (&charN, &nm2, &n, &scalar, grid->get_data (2) + 1, &n, data_in, &ione, &dpone, data_out + 1, &ione);

			data_out [n - 1] += scalar * ddot_ (&n, grid->get_data (2) + n - 1, &n, data_in, &ione);

			TRACE (logger, "Operation complete.");
		}

		implicit_diffusion::implicit_diffusion (double i_coeff, double& i_timestep, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_matrix, int *i_flags_ptr, int i_logger) : bases::implicit_plan (i_n, i_grid, i_matrix, i_flags_ptr, i_logger), timestep (i_timestep) {
			coeff = i_coeff;
			timestep = i_timestep;
			n = i_n;
			grid = i_grid;
			matrix = i_matrix;
		}

		void implicit_diffusion::execute () {
			bases::implicit_plan::execute ();
						
			TRACE (logger, "Operating...");
					
			double scalar = coeff * timestep;
			// This is the main loop for setting up the diffusion equation in Chebyshev space
			for (int i = 1; i < n - 1; ++i) {
			   	daxpy_ (&n, &scalar, grid->get_data (2) + i, &n, matrix + i, &n);
			}
	
			TRACE (logger, "Operation complete.");
		}
	} /* chebyshev */
} /* one_d */

