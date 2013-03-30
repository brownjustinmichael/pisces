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
	cheb_1D::cheb_1D (double i_coeff, int i_n, double *i_data_in, double *i_data_out) {
		int i, j;
		coeff = i_coeff;
		n = i_n;
		data_in = i_data_in;
		if (i_data_out == NULL) {
			data_out = i_data_in;	
		} else {
			data_out = i_data_out;
		}
		
		LOG4CXX_DEBUG (config::logger, "cheb_1D: [cheb_1d] Input pointer: " << data_in << ", Output pointer: " << data_out);
		
		even_diffusion_matrix.resize (i_n * (i_n + 6) / 8);
		odd_diffusion_matrix.resize (i_n * (i_n + 6) / 8);
		
		LOG4CXX_TRACE (config::logger, "cheb_1D: [cheb_1d] Instantiating...");
		
	   	for (j = 0; j < i_n / 2; ++j) {
			even_diffusion_matrix [j*(j+1)/2] = -4*j*j*j*coeff;
		   	odd_diffusion_matrix [j*(j+1)/2] = (2*j+1)*(1-pow(2*j+1,2))*coeff;
		   	for (i=1;i<=j;++i) {
				even_diffusion_matrix [i+j*(j+1)/2] = 8*j*(i*i-j*j)*coeff;
			   	odd_diffusion_matrix [i+j*(j+1)/2] = (2*j+1)*(pow(2*i+1,2)-pow(2*j+1,2))*coeff;
			}
		}
		
		LOG4CXX_TRACE (config::logger, "cheb_1D: [cheb_1d] Instantiation complete.");
	}

	void cheb_1D::execute (double timestep) {
		int i;
	    int ione=1, itwo=2, nhalf = n / 2;
	    char charN = 'N', charU = 'U';
	    double dpone = 1.e0, dmone = -1.e0;
	
		LOG4CXX_TRACE (config::logger, "cheb_1D: [operate] Operating...");
		
		for (i = 0; i < n / 2; ++i) {
			even_diffusion_matrix [i*(i+3)/2] = 1.0 / timestep * coeff;
			odd_diffusion_matrix [i*(i+3)/2] = 1.0 / timestep * coeff;
		}
	
		if (data_in != data_out) {
			dcopy_ (&n, data_in, &ione, data_out, &ione);
		}
	
	    dtpsv_(&charU, &charN, &charN, &nhalf, &(even_diffusion_matrix [0]), data_out, &itwo);
	    dtpsv_(&charU, &charN, &charN, &nhalf, &(odd_diffusion_matrix [0]), &(data_out [1]), &itwo);
		
		for (i = 0; i < n; ++i) {
			data_out [i] /= timestep;
		}
		
		LOG4CXX_TRACE (config::logger, "cheb_1D: [operate] Operation complete.");
	}
	
} /* diffusion */
