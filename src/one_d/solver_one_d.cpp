/*!***********************************************************************
 * \file one_d/solver.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-13.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "../config.hpp"
#include "solver_one_d.hpp"

namespace one_d
{
	void lapack_solver::i_factorize () {
		int info;
		dgetrf_ (&n, &n, matrix, &n, &ipiv [0], &info);
		
		if (info != 0) {
			ERROR (logger, "Unable to invert matrix");
			throw 0; // For now, kill the program. 
			/*
				TODO Replace this with a more useful exception that can be handled
			*/
		}
	}

	void lapack_solver::i_solve () {
		int ione = 1, info;
		double dpone = 1.0;
		char charN = 'N';
		
		TRACE (logger, "Solving...")
		
		if (data_out == rhs) {
			daxpy_ (&n, &dpone, data_in, &ione, data_out, &ione);
		} else {
			if (data_in != data_out) {
				dcopy_ (&n, data_in, &ione, data_out, &ione);
			}
			daxpy_ (&n, &dpone, rhs, &ione, data_out, &ione);
		}
		
		dgetrs_ (&charN, &n, &ione, &matrix [0], &n, &ipiv [0], data_out, &n, &info);
		
		if (info != 0) {
			ERROR (logger, "Unable to solve factorized matrix equation");
			throw 0; // For now, kill the program. 
			/*
				TODO Replace this with a more useful exception that can be handled
			*/
		}
		
		TRACE (logger, "Solve complete.")
	}
} /* one_d */