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
	collocation_chebyshev_implicit_1D::collocation_chebyshev_implicit_1D (double i_coeff, int i_n, std::shared_ptr<collocation::chebyshev_grid> i_cheb, double *i_matrix, int i_flags) {
		coeff = i_coeff;
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

	void collocation_chebyshev_implicit_1D::execute (double timestep) {
	    int ione = 1, info, i, j, start, end;
	    char charN = 'N';
	    double dpone = 1.e0, dzero = 0.0;
			
		TRACE ("Operating...");

		start = 0;
		end = n;
		if (flags & boundary::fixed_upper) {
			// This equation fixes the value at the top boundary
			start = 1;
		}
	
		if (flags & boundary::fixed_lower) {
			// This equation fixes the value at the bottom boundary
			end = n - 1;
		}
	
		// This is the main loop for setting up the diffusion equation in Chebyshev space
		for (i = start; i < end; ++i) {
		   	for (j = 0; j < n; ++j) {
				matrix [i + j * n] += coeff * timestep * cheb->index (2, j, i);
			}
		}

		previous_timestep = timestep;
		
		TRACE ("Operation complete.");

	}
	
	collocation_chebyshev_explicit_1D::collocation_chebyshev_explicit_1D (double i_coeff, int i_n, std::shared_ptr<collocation::chebyshev_grid> i_cheb, double *i_data_in, double *i_data_out, int i_flags) {
		coeff = i_coeff;
		n = i_n;
		data_in = i_data_in;
		data_out = i_data_out;
		
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

	void collocation_chebyshev_explicit_1D::execute (double timestep) {
	    int ione = 1;
	    char charN = 'N';
	    double dpone = 1.0;
		
		TRACE ("Operating...");
		
		double scalar = coeff * timestep;
		
		// Set up and evaluate the explicit part of the diffusion equation
		dgemv_ (&charN, &n, &n, &dpone, cheb->get_data (2), &n, &data_in [0], &ione, &dpone, &data_out [0], &ione);
		
		if (flags & boundary::fixed_upper) {
			data_out [0] = 0.0;
		}
		if (flags & boundary::fixed_lower) {
			data_out [n - 1] = 0.0;			
		}
		
		dscal_ (&n, &scalar, &data_out [0], &ione);
		
		TRACE ("Operation complete.");

	}
} /* diffusion */
