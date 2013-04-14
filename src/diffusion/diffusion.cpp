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
	collocation_chebyshev_implicit_1D::collocation_chebyshev_implicit_1D (double i_coeff, double i_alpha, int i_n, std::shared_ptr<collocation::chebyshev_grid> i_cheb, double *i_data_in, double *i_rhs, double *i_data_out, int i_flags) {
		coeff = i_coeff;
		alpha = i_alpha;
		previous_timestep = 0.0;
		n = i_n;
		
		assert (i_data_in != i_rhs);
		
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
			cheb.reset (new collocation::chebyshev_grid (i_n, i_n, sqrt (2.0 / (i_n - 1))));
		} else {
			TRACE ("Existing collocation grid, sharing pointer...")
			cheb = i_cheb;
		}
		
		DEBUG ("Input pointer: " << data_in << ", Output pointer: " << data_out);
		
		diffusion_matrix.resize (i_n * i_n);
		ipiv.resize (i_n * i_n);
		pre_matrix.resize (i_n * i_n);
		temp.resize (i_n);
		
		TRACE ("Instantiation complete.");
	}

	void collocation_chebyshev_implicit_1D::execute (double timestep) {
	    int ione = 1, info, i, j;
	    char charN = 'N';
	    double dpone = 1.e0, dzero = 0.0;
			
		TRACE ("Operating...");
				
		// if (data_out == rhs) {
		// 	dscal_ (&n, &timestep, &rhs [0], &ione);
		// 	daxpy_ (&n, &dpone, &data_in [0], &ione, &data_out [0], &ione);
		// } else {
		// 	if (data_in != data_out) {
		// 		dcopy_ (&n, &data_in [0], &ione, &data_out [0], &ione);
		// 	}
		// 	daxpy_ (&n, &timestep, &rhs [0], &ione, &data_out [0], &ione);
		// }
		// 		
		
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
	
	void collocation_chebyshev_implicit_1D::matrix (double alpha_scalar, double *matrix) {
		int i, j, start, end;
		
		start = 0;
		end = n;
		if (flags & boundary::fixed_upper) {
			// This equation fixes the value at the top boundary
			start = 1;
			for (j = 0; j < n; ++j) {
				matrix [j * n] = cheb->index (0, j, 0);
			}
		}
		
		if (flags & boundary::fixed_lower) {
			// This equation fixes the value at the bottom boundary
			end = n - 1;
			for (j = 0; j < n; ++j) {
				matrix [(n - 1) + j * n] = cheb->index (0, j, n - 1);
			}
		}
		
		// This is the main loop for setting up the diffusion equation in Chebyshev space
		for (i = start; i < end; ++i) {
		   	for (j = 0; j < n; ++j) {
				matrix [i + j * n] = (cheb->index (0, j, i) + alpha_scalar * cheb->index (2, j, i));
			}
		}
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
	    int i, ione = 1;
	    char charN = 'N';
	    double dpone = 1.0, dzero = 0.0;
		
		TRACE ("Operating...");
		
		double scalar = coeff * timestep;
		
		// Set up and evaluate the explicit part of the diffusion equation
		dgemv_ (&charN, &n, &n, &dpone, cheb->get_data (2), &n, &data_in [0], &ione, &dzero, &data_out [0], &ione);
		data_out [0] = 0.0;
		data_out [n - 1] = 0.0;
		dscal_ (&n, &scalar, &data_out [0], &ione);
		
		TRACE ("Operation complete.");

	}
} /* diffusion */
