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
		if (i_data_out == NULL) {
			data_out = i_data_in;
		} else {
			data_out = i_data_out;
		}

		LOG4CXX_TRACE (config::logger, "Instantiating...");

		cheb = new collocation::cheb_grid (i_n, i_n);
		
		LOG4CXX_DEBUG (config::logger, "Input pointer: " << data_in << ", Output pointer: " << data_out);

		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) {
				LOG4CXX_DEBUG (config::logger, "cheb [" << i << ", " << j << "] = " << cheb->index (0, i, j));
			}
		}
		
		diffusion_matrix.resize (i_n * i_n);
		
		LOG4CXX_TRACE (config::logger, "Instantiation complete.");
	}

	void cheb_1D::execute (double timestep) {
		int i, j;
	    int ione=1, itwo=2, nhalf = n / 2;
	    char charN = 'N', charU = 'U';
	    double dpone = 1.e0, dmone = -1.e0;
	
		LOG4CXX_TRACE (config::logger, "Operating...");
		
		double scalar = std::sqrt (2.0 / (n - 1.0));
		double scalar_half = scalar / 2.0;
		double alpha_scalar = alpha * timestep * coeff;
	   	for (i = 0; i < n; ++i) {
			diffusion_matrix [i * n] = scalar_half * (cheb->index (0, i, 0) - alpha_scalar * cheb->index (2, i, 0));
		   	for (j = 1; j < n - 1; ++j) {
				diffusion_matrix [i * n + j] = scalar * (cheb->index (0, i, j) - alpha_scalar * cheb->index (2, i, j));
			}
			diffusion_matrix [i * n + n - 1] = scalar_half * (cheb->index (0, i, n - 1) - alpha_scalar * cheb->index (2, i, n - 1));
		}
	
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
