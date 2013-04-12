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
	collocation_chebyshev_1D::collocation_chebyshev_1D (double i_coeff, double i_alpha, int i_n, std::shared_ptr<collocation::chebyshev_grid> i_cheb, double *i_data_in, double *i_rhs, double *i_data_out, int i_flags) {
		coeff = i_coeff;
		alpha = i_alpha;
		previous_timestep = 0.0;
		n = i_n;
		data_in = i_data_in;
		rhs = i_rhs;
		if (i_data_out == NULL) {
			data_out = i_data_in;
		} else {
			data_out = i_data_out;
		}
		
		flags = i_flags;

		TRACE ("Instantiating...");

		if (!i_cheb) {
			TRACE ("No collocation grid yet, constructing...")
			cheb.reset (new collocation::chebyshev_grid (i_n, i_n));
		} else {
			TRACE ("Existing collocation grid, sharing pointer...")
			cheb = i_cheb;
		}
		
		DEBUG ("Input pointer: " << data_in << ", Output pointer: " << data_out);
		
		diffusion_matrix.resize (i_n * i_n);
		ipiv.resize (i_n * i_n);
		pre_matrix.resize (i_n * i_n);
		
		TRACE ("Instantiation complete.");
	}

	void collocation_chebyshev_1D::execute (double timestep) {
	    int ione = 1, info;
	    char charN = 'N';
	    double dpone = 1.e0, dzero = 0.0;
			
		TRACE ("Operating...");
		
		// Set up and evaluate the explicit part of the diffusion equation
		if (timestep != previous_timestep) {
			matrix ((1.0 - alpha) * timestep * coeff, &pre_matrix [0]);
		}
		dgemv_ (&charN, &n, &n, &dpone, &pre_matrix [0], &n, &data_in [0], &ione, &dpone, &rhs [0], &ione);

		if (rhs != data_out) {
			dcopy_ (&n, &rhs [0], &ione, &data_out [0], &ione);
		}

		// Set up and evaluate the implicit part of the diffusion equation
		if (timestep != previous_timestep) {
			matrix (- alpha * timestep * coeff, &diffusion_matrix [0]);
			dgetrf_ (&n, &n, &diffusion_matrix [0], &n, &ipiv [0], &info);			
		}
		dgetrs_ (&charN, &n, &ione, &diffusion_matrix [0], &n, &ipiv [0], &data_out [0], &n, &info);
		
		if (info != 0) {
			ERROR ("Unable to invert matrix");
			// throw 0; // For now, kill the program. 
			/*
				TODO Replace this with a more useful exception that can be handled
			*/
		}
		
		previous_timestep = timestep;
		
		TRACE ("Operation complete.");

	}
	
	void collocation_chebyshev_1D::matrix (double alpha_scalar, double *matrix) {
		int i, j, start, end;
		double scalar = std::sqrt (2.0 / (n - 1.0));
		double scalar_half = scalar / 2.0;
		
		start = 0;
		end = n;
		if (flags & boundary::fixed_upper) {
			// This equation fixes the value at the top boundary
			start = 1;
			matrix [0] = scalar_half * cheb->index (0, 0, 0);
			for (j = 1; j < n - 1; ++j) {
				matrix [j * n] = scalar * cheb->index (0, j, 0);
			}
			matrix [(n - 1) * n] = scalar_half * cheb->index (0, n - 1, 0);
		}
		
		if (flags & boundary::fixed_lower) {
			// This equation fixes the value at the bottom boundary
			end = n - 1;
			matrix [(n - 1)] = scalar_half * cheb->index (0, 0, n - 1);
			for (j = 1; j < n - 1; ++j) {
				matrix [(n - 1) + j * n] = scalar * cheb->index (0, j, n - 1);
			}
			matrix [(n - 1) + (n - 1) * n] = scalar_half * cheb->index (0, n - 1, n - 1);
		}
		
		// This is the main loop for setting up the diffusion equation in Chebyshev space
		for (i = start; i < end; ++i) {
			matrix [i] = scalar_half * (cheb->index (0, 0, i) + alpha_scalar * cheb->index (2, 0, i));
		   	for (j = 1; j < n - 1; ++j) {
				matrix [i + j * n] = scalar * (cheb->index (0, j, i) + alpha_scalar * cheb->index (2, j, i));
			}
			matrix [i + (n - 1) * n] = scalar_half * (cheb->index (0, n - 1, i) + alpha_scalar * cheb->index (2, n - 1, i));
		}
	}
} /* diffusion */
