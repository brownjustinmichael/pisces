/*!***********************************************************************
 * \file solver.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-13.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "../config.hpp"
#include "solver.hpp"

namespace solver
{
	void lapack_solver::solve (int *flags) {
		int ione = 1, info;
		double dpone = 1.0;
		char charN = 'N';
		
		if (data_out == rhs) {
			daxpy_ (&n, &dpone, &data_in [0], &ione, &data_out [0], &ione);
		} else {
			if (data_in != data_out) {
				dcopy_ (&n, &data_in [0], &ione, &data_out [0], &ione);
			}
			daxpy_ (&n, &dpone, &rhs [0], &ione, &data_out [0], &ione);
		}
		
		// Set up and evaluate the implicit part of the diffusion equation
		if (! (*flags & factorized)) {
			dgetrf_ (&n, &n, &matrix [0], &n, &ipiv [0], &info);
			*flags |= factorized; 
		}
				
		dgetrs_ (&charN, &n, &ione, &matrix [0], &n, &ipiv [0], &data_out [0], &n, &info);
		
		if (info != 0) {
			ERROR ("Unable to invert matrix");
			throw 0; // For now, kill the program. 
			/*
				TODO Replace this with a more useful exception that can be handled
			*/
		}
	}
} /* solver */