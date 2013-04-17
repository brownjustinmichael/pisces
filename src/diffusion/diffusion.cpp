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
	namespace implicit_methods
	{
		collocation_chebyshev_1D::collocation_chebyshev_1D (double i_coeff, double *i_timestep_ptr, int i_n, std::shared_ptr<collocation::chebyshev_grid> i_cheb, double *i_matrix, int *i_flags) {
			coeff = i_coeff;
			timestep_ptr = i_timestep_ptr;
			previous_timestep = 0.0;
			n = i_n;
			matrix = i_matrix;
		
			flags = i_flags;

			TRACE ("Instantiating...");

			if (!i_cheb) {
				TRACE ("No collocation grid yet, constructing...")
				cheb.reset (new collocation::chebyshev_grid (i_n, i_n, sqrt (2.0 / (i_n - 1))));
			} else {
				TRACE ("Existing collocation grid, sharing pointer...")
				cheb = i_cheb;
			}
				
			TRACE ("Instantiation complete.");
		}

		void collocation_chebyshev_1D::execute () {
		    int ione = 1, i, start, len;
			double scalar;
			
			TRACE ("Operating...");

			if (!flags) {
				start = 0;
			} else if (*flags & boundary::fixed_upper) {
				// This equation fixes the value at the top boundary
				start = 1;
			} else {
				start = 0;
			}
	
			if (!flags) {
				len = n - start;
			} else if (*flags & boundary::fixed_lower) {
				// This equation fixes the value at the bottom boundary
				len = n - 1 - start;
			} else {
				len = n - start;
			}
		
			scalar = coeff * *timestep_ptr;
	
			// This is the main loop for setting up the diffusion equation in Chebyshev space
			for (i = 0; i < n; ++i) {
			   	daxpy_ (&len, &scalar, cheb->get_data (2) + start + i * n, &ione, &matrix [start + i * n], &ione);
			}
		
			TRACE ("Operation complete.");
		}
	} /* implicit */

	namespace explicit_methods
	{
		collocation_chebyshev_1D::collocation_chebyshev_1D (double i_coeff, double *i_timestep_ptr, int i_n, std::shared_ptr<collocation::chebyshev_grid> i_cheb, double *i_data_in, double *i_data_out, int *i_flags) {
			coeff = i_coeff;
			timestep_ptr = i_timestep_ptr;
			n = i_n;
			data_in = i_data_in;
			data_out = i_data_out;
		
			if (!i_flags) {
				default_flags = 0x00;
				flags = &default_flags;
			} else {
				flags = i_flags;
			}

			TRACE ("Instantiating...");

			if (!i_cheb) {
				TRACE ("No collocation grid yet, constructing...")
				cheb.reset (new collocation::chebyshev_grid (i_n, i_n, sqrt (2.0 / (i_n - 1))));
			} else {
				TRACE ("Existing collocation grid, sharing pointer...")
				cheb = i_cheb;
			}

			TRACE ("Instantiation complete.");
		}

		void collocation_chebyshev_1D::execute () {
		    int ione = 1, len, start;
		    char charN = 'N';
		    double dpone = 1.0, scalar;
		
			TRACE ("Operating...");
		
			scalar = coeff * *timestep_ptr;
			
			for (int i = 0; i < n; ++i) {
				DEBUG ("data [" << i << "] = " << data_out [i]);
			}		
			// Set up and evaluate the explicit part of the diffusion equation
			dgemv_ (&charN, &n, &n, &scalar, cheb->get_data (2), &n, data_in, &ione, &dpone, data_out, &ione);
			
			for (int i = 0; i < n; ++i) {
				DEBUG ("data [" << i << "] = " << data_out [i]);
			}

			TRACE ("Operation complete.");
		}
	} /* explicit */
} /* diffusion */
