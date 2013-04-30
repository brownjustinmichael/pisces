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
#include "../utils/chebyshev.hpp"
#include "../utils/utils.hpp"

namespace one_d
{
	namespace chebyshev
	{
		void explicit_diffusion::init (double i_coeff, std::shared_ptr<bases::collocation_grid> i_grid, double* i_position) {
			coeff = i_coeff;
			grid = i_grid;
			position = i_position;
		}

		void explicit_diffusion::execute () {
			double scalar = coeff * timestep;
			double dy0, dyn;
			
			bases::explicit_plan::execute ();
		
			TRACE (logger, "Operating...");
		
			// dy0 = scalar * utils::dot (n, grid->get_data (1) + 1, data_in, n);
			// data_out [0] += dy0;
			// dyn = scalar * utils::dot (n, grid->get_data (0) + n - 1, data_in, n);
			// data_out [n - 1] += dyn;
			data_out [0] += scalar * utils::dot (n, grid->get_data (2), data_in, n);
			
			// Set up and evaluate the explicit part of the diffusion equation
			utils::matrix_vector_multiply (n - 2, n, scalar, grid->get_data (2) + 1, data_in, 1.0, data_out + 1, n);

			data_out [n - 1] += scalar * utils::dot (n, grid->get_data (2) + n - 1, data_in, n);

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
			   	utils::add_scaled (n, scalar, grid->get_data (2) + i, matrix + i, n, n);
			}
	
			TRACE (logger, "Operation complete.");
		}
	} /* chebyshev */
} /* one_d */

