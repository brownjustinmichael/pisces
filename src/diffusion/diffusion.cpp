// 
//! \file diffusion.cpp
//  src
//  
//  Created by Justin Brown on 2013-03-22.
//  Copyright 2013 Justin Brown. All rights reserved.
// 

#include <iostream>
#include <cmath>
#include <cassert>
#include "../config.hpp"
#include "diffusion.hpp"
#include "../collocation/collocation.hpp"
#include "../boundary/boundary.hpp"

namespace diffusion
{
	cheb_1D::cheb_1D (double i_coeff, double i_alpha, int i_n, double *i_data_in, double *i_data_out, int i_flags) {
		int i, j;
		coeff = i_coeff;
		alpha = i_alpha;
		n = i_n;
		data_in = i_data_in;
		data_in = i_data_in;
		if (i_data_out == NULL) {
			data_out = i_data_in;
		} else {
			data_out = i_data_out;
		}
		flags = i_flags;

		TRACE ("Instantiating...");

		cheb = new collocation::cheb_grid (i_n, i_n);
		
		DEBUG ("Input pointer: " << data_in << ", Output pointer: " << data_out);
		
		diffusion_matrix.resize (i_n * i_n);
		ipiv.resize (i_n * i_n);
		pre_matrix.resize (i_n * i_n);
		temp.resize (i_n);
		
		TRACE ("Instantiation complete.");
	}

	void cheb_1D::execute (double timestep) {
		int i, j;
	    int ione=1, itwo=2, info;
	    char charN = 'N', charU = 'U';
	    double dpone = 1.e0, dmone = -1.e0, dzero = 0.0;
		double d2sum;
	
		TRACE ("Operating...");
		
		dcopy_ (&n, &data_in [0], &ione, &temp [0], &ione);
		
		// Set up and evaluate the explicit part of the diffusion equation
		matrix ((1.0 - alpha) * timestep * coeff, &pre_matrix [0]);
		dgemv_ (&charN, &n, &n, &dpone, &pre_matrix [0], &n, &temp [0], &ione, &dzero, &data_out [0], &ione);

		// Set up and evaluate the implicit part of the diffusion equation
		matrix (- alpha * timestep * coeff, &diffusion_matrix [0]);
		
		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) {
				DEBUG ("matrix [" << i << ", " << j << "] = " << diffusion_matrix [i + j * n]);
			}
		}
		for (i = 0; i < n; ++i) {
			DEBUG ("in [" << i << "] = " << data_out [i]);
		}
		dgesv_ (&n, &ione, &diffusion_matrix [0], &n, &ipiv [0], &data_out [0], &n, &info);
		for (i = 0; i < n; ++i) {
			DEBUG ("out [" << i << "] = " << data_out [i]);
		}	
		TRACE ("Operation complete.");
	}
	
	void cheb_1D::matrix (double alpha_scalar, double *matrix) {
		int i, j;
		double scalar = std::sqrt (2.0 / (n - 1.0));
		double scalar_half = scalar / 2.0;
		
		// This equation fixes the value at the top boundary
		matrix [0] = scalar_half * cheb->index (0, 0, 0);
		for (j = 1; j < n - 1; ++j) {
			matrix [j * n] = scalar * cheb->index (0, j, 0);
		}
		matrix [(n - 1) * n] = scalar_half * cheb->index (0, n - 1, 0);
		
		// This is the main loop for setting up the diffusion equation in Chebyshev space
		for (i = 1; i < n - 1; ++i) {
			matrix [i] = scalar_half * (cheb->index (0, 0, i) + alpha_scalar * cheb->index (2, 0, i));
		   	for (j = 1; j < n - 1; ++j) {
				matrix [i + j * n] = scalar * (cheb->index (0, j, i) + alpha_scalar * cheb->index (2, j, i));
			}
			matrix [i + (n - 1) * n] = scalar_half * (cheb->index (0, n - 1, i) + alpha_scalar * cheb->index (2, n - 1, i));
		}
		
		// This equation fixes the value at the bottom boundary
		matrix [(n - 1)] = scalar_half * cheb->index (0, 0, n - 1);
		for (j = 1; j < n - 1; ++j) {
			matrix [(n - 1) + j * n] = scalar * cheb->index (0, j, n - 1);
		}
		matrix [(n - 1) + (n - 1) * n] = scalar_half * cheb->index (0, n - 1, n - 1);
	}
} /* diffusion */
