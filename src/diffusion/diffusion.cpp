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

namespace diffusion
{
	cheb_1D::cheb_1D (double i_coeff, double i_alpha, int i_n, double *i_data_in, double *i_data_out) {
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

		LOG4CXX_TRACE (config::logger, "Instantiating...");

		cheb = new collocation::cheb_grid (i_n, i_n);
		
		LOG4CXX_DEBUG (config::logger, "Input pointer: " << data_in << ", Output pointer: " << data_out);
		
		diffusion_matrix.resize (i_n * i_n);
		ipiv.resize (i_n * i_n);
		pre_matrix.resize (i_n * i_n);
		temp.resize (i_n);
		
		LOG4CXX_TRACE (config::logger, "Instantiation complete.");
	}

	void cheb_1D::execute (double timestep) {
		int i, j;
	    int ione=1, itwo=2, info;
	    char charN = 'N', charU = 'U';
	    double dpone = 1.e0, dmone = -1.e0, dzero = 0.0;
		double d2sum;
	
		LOG4CXX_TRACE (config::logger, "Operating...");
		
		dcopy_ (&n, &data_in [0], &ione, &temp [0], &ione);
		
		double scalar = std::sqrt (2.0 / (n - 1.0));
		double scalar_half = scalar / 2.0;
		double alpha_scalar = (1.0 - alpha) * timestep * coeff;
		
		pre_matrix [0] = scalar_half * cheb->index (0, 0, 0);
		for (j = 1; j < n - 1; ++j) {
			pre_matrix [j * n] = scalar * cheb->index (0, j, 0);
		}
		pre_matrix [(n - 1) * n] = scalar_half * cheb->index (0, n - 1, 0);
		for (i = 1; i < n - 1; ++i) {
			pre_matrix [i] = scalar_half * (cheb->index (0, 0, i) + alpha_scalar * cheb->index (2, 0, i));
		   	for (j = 1; j < n - 1; ++j) {
				pre_matrix [i + j * n] = scalar * (cheb->index (0, j, i) + alpha_scalar * cheb->index (2, j, i));
			}
			pre_matrix [i + (n - 1) * n] = scalar_half * (cheb->index (0, n - 1, i) + alpha_scalar * cheb->index (2, n - 1, i));
		}
		pre_matrix [(n - 1)] = scalar_half * cheb->index (0, 0, n - 1);
		for (j = 1; j < n - 1; ++j) {
			pre_matrix [(n - 1) + j * n] = scalar * cheb->index (0, j, n - 1);
		}
		pre_matrix [(n - 1) + (n - 1) * n] = scalar_half * cheb->index (0, n - 1, n - 1);
			
		alpha_scalar = alpha * timestep * coeff;
		
		diffusion_matrix [0] = scalar_half * cheb->index (0, 0, 0);
		for (j = 1; j < n - 1; ++j) {
			diffusion_matrix [j * n] = scalar * cheb->index (0, j, 0);
		}
		diffusion_matrix [(n - 1) * n] = scalar_half * cheb->index (0, n - 1, 0);
		for (i = 1; i < n - 1; ++i) {
			diffusion_matrix [i] = scalar_half * (cheb->index (0, 0, i) - alpha_scalar * cheb->index (2, 0, i));
		   	for (j = 1; j < n - 1; ++j) {
				diffusion_matrix [i + j * n] = scalar * (cheb->index (0, j, i) - alpha_scalar * cheb->index (2, j, i));
			}
			diffusion_matrix [i + (n - 1) * n] = scalar_half * (cheb->index (0, n - 1, i) - alpha_scalar * cheb->index (2, n - 1, i));
		}
		diffusion_matrix [(n - 1)] = scalar_half * cheb->index (0, 0, n - 1);
		for (j = 1; j < n - 1; ++j) {
			diffusion_matrix [(n - 1) + j * n] = scalar * cheb->index (0, j, n - 1);
		}
		diffusion_matrix [(n - 1) + (n - 1) * n] = scalar_half * cheb->index (0, n - 1, n - 1);
		
		dgemv_ (&charN, &n, &n, &dpone, &pre_matrix [0], &n, &temp [0], &ione, &dzero, &data_out [0], &ione);

		dgesv_ (&n, &ione, &diffusion_matrix [0], &n, &ipiv [0], &data_out [0], &n, &info);
	
		LOG4CXX_TRACE (config::logger, "Operation complete.");
	}
	
	fixed_angle_1D::fixed_angle_1D (double i_coeff, int i_n, double *i_data_in, double *i_data_out) {
		int i;
		coeff = i_coeff;
		n = i_n;
		data_in = i_data_in;
		if (i_data_out == NULL) {
			data_out = i_data_in;	
		} else {
			data_out = i_data_out;
		}
		
		diagonal.resize (i_n);
		diagonal_right.resize (i_n - 1);
		diagonal_left.resize (i_n - 1);
	}
	
	void fixed_angle_1D::execute (double timestep) {
		int i;
		double pi = std::acos (-1.0);
		double angle = pi / n;
		double scalar = coeff * timestep * n / 4.0 / pi;
		double scalartwo = scalar * 2.0;
		double sine;
		int ione = 1;
		int info;
		diagonal [0] = 1.0;
		diagonal_right [0] = 0.0;
		for (i = 1; i < n - 1; ++i) {
			sine = std::sin (angle * i);
			diagonal [i] = scalartwo / sine / sine + 1.0;
			diagonal_right [i] = - scalar / sine / sine;
			diagonal_left [i - 1] = diagonal_right [i];
		}
		diagonal [n - 1] = 1.0;
		diagonal_left [n - 2] = 0.0;
		
		dgtsv_ (&n, &ione, &diagonal_left [0], &diagonal [0], &diagonal_right [0], data_out, &n, &info);
	}
} /* diffusion */
