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

namespace diffusion
{
	plan::plan (operation *i_operation_ptr, int i_n_data_ptrs, double **i_data_ptrs) {
		int i;
	
		// Check that the number of data pointers for the operation equals the number passed
		assert (i_operation_ptr->n_data_ptrs == i_n_data_ptrs);

		operation_ptr = i_operation_ptr;

		LOG4CXX_TRACE (config::logger, "plan: Instantiating...");

		for (i = 0; i < i_n_data_ptrs; ++i) {
			data_ptrs.push_back (i_data_ptrs [i]);
		}
	
		LOG4CXX_TRACE (config::logger, "plan: Instantiation complete.");
	}

	cheb_1D::cheb_1D (int i_n, double i_coeff) : operation (i_n, 2) {
		int i, j;
		coeff = i_coeff;
		even_diffusion_matrix.resize (i_n / 2*(i_n / 2+3)/2);
		odd_diffusion_matrix.resize (i_n / 2*(i_n / 2+3)/2);
		
		LOG4CXX_TRACE (config::logger, "cheb_1D: Instantiating...");
		
	   	for (j = 0; j < i_n / 2; ++j) {
			even_diffusion_matrix [j*(j+1)/2] = -4*j*j*j*coeff;
		   	odd_diffusion_matrix [j*(j+1)/2] = (2*j+1)*(1-pow(2*j+1,2))*coeff;
		   	for (i=1;i<=j;++i) {
				even_diffusion_matrix [i+j*(j+1)/2] = 8*j*(i*i-j*j)*coeff;
			   	odd_diffusion_matrix [i+j*(j+1)/2] = (2*j+1)*(pow(2*i+1,2)-pow(2*j+1,2))*coeff;
			}
		}
		
		LOG4CXX_TRACE (config::logger, "cheb_1D: Instantiation complete.");
	}

	void cheb_1D::operate (double timestep, double **data_ptrs) {
		int i;
	    int ione=1, itwo=2, nhalf = n / 2;
	    char charN = 'N', charU = 'U';
	    double dpone = 1.e0, dmone = -1.e0;
	
		LOG4CXX_TRACE (config::logger, "cheb_1D: [operate] Operating...");
		
		for (i = 0; i < n / 2; ++i) {
			even_diffusion_matrix [i*(i+3)/2] = 1.0 / timestep * coeff;
			odd_diffusion_matrix [i*(i+3)/2] = 1.0 / timestep * coeff;
		}
	
		if (data_ptrs [0] != data_ptrs [1]) {
			dcopy_ (&n, data_ptrs [0], &ione, data_ptrs [1], &ione);
		}
	
	    dtpsv_(&charU, &charN, &charN, &nhalf, &(even_diffusion_matrix [0]), data_ptrs [1], &itwo);
	    dtpsv_(&charU, &charN, &charN, &nhalf, &(odd_diffusion_matrix [0]), &(data_ptrs [1] [1]), &itwo);
		
		for (i = 0; i < n; ++i) {
			data_ptrs [1] [i] /= timestep;
		}
		
		LOG4CXX_TRACE (config::logger, "cheb_1D: [operate] Operation complete.");
	}
	
} /* diffusion */
