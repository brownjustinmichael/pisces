/*!***********************************************************************
 * \file diffusion.cpp
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
#include "../boundary/boundary.hpp"

namespace diffusion
{
	namespace explicit_methods
	{
		collocation_chebyshev_1D::collocation_chebyshev_1D (double i_coeff, double *i_timestep_ptr, int i_n, std::shared_ptr<collocation::chebyshev_grid> i_grid, double *i_data_in, double *i_data_out, int *i_flags) {
			coeff = i_coeff;
			timestep_ptr = i_timestep_ptr;
			n = i_n;
			data_in = i_data_in;
			data_out = i_data_out;
			
			TRACE ("Instantiating...");
		
			if (!i_flags) {
				default_flags = 0x00;
				flags = &default_flags;
			} else {
				flags = i_flags;
			}

			if (!i_grid) {
				TRACE ("No collocation grid yet, constructing...")
				grid.reset (new collocation::chebyshev_grid (i_n, i_n, sqrt (2.0 / (i_n - 1))));
			} else {
				TRACE ("Existing collocation grid, sharing pointer...")
				grid = i_grid;
			}

			TRACE ("Instantiation complete.");
		}

		void collocation_chebyshev_1D::execute () {
		    int ione = 1;
		    char charN = 'N';
		    double dpone = 1.0;
		
			TRACE ("Operating...");
		
			double scalar = coeff * *timestep_ptr;
			// Set up and evaluate the explicit part of the diffusion equation
			dgemv_ (&charN, &n, &n, &scalar, grid->get_data (2), &n, data_in, &ione, &dpone, data_out, &ione);

			TRACE ("Operation complete.");
		}
	} /* explicit */
	
	namespace implicit_methods
	{
		collocation_chebyshev_1D::collocation_chebyshev_1D (double i_coeff, double i_alpha_0, double i_alpha_n, double *i_timestep_ptr, int i_n, std::shared_ptr<collocation::chebyshev_grid> i_grid, double *i_matrix, int *i_flags) {
			coeff = i_coeff;
			alpha_0 = i_alpha_0;
			alpha_n = i_alpha_n;
			timestep_ptr = i_timestep_ptr;
			n = i_n;
			matrix = i_matrix;
		
			flags = i_flags;

			TRACE ("Instantiating...");

			if (!i_grid) {
				TRACE ("No collocation grid yet, constructing...")
				grid.reset (new collocation::chebyshev_grid (i_n, i_n, sqrt (2.0 / (i_n - 1))));
			} else {
				TRACE ("Existing collocation grid, sharing pointer...")
				grid = i_grid;
			}
				
			TRACE ("Instantiation complete.");
		}

		void collocation_chebyshev_1D::execute () {			
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
	} /* implicit */
} /* diffusion */
