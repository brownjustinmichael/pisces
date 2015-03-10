/*!**********************************************************************
 * \file tridiagonal.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef TRIDIAGONAL_CPP_53CA7D2B
#define TRIDIAGONAL_CPP_53CA7D2B

#include <vector>
#ifdef _MPI
#include <mpi.h>
#endif
#include "linalg/utils.hpp"
#include "linalg/linalg.hpp"
#include "tridiagonal.hpp"

namespace linalg
{
	namespace block
	{
		void tridiag_factorize (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, int* ipiv, double *x, int *xipiv, int *info, int nrhs, int lda) {
			int ntop, nbot;
			int ldx = 2;
			if (id == 0) {
				ntop = 0;
			} else {
				ntop = 1;
			}
			if (id == np - 1) {
				nbot = 0;
			} else {
				nbot = 1;
			}
		
			if (lda == -1) {
				lda = n + ntop + nbot;
			}
		
			int ldxx = 2 * nrhs;
			double *bufferl = x, *bufferr = bufferl + n * nrhs, *xsub = bufferr + n * nrhs, *xdiag = xsub + ldxx, *xsup = xdiag + ldxx, *xsupsup = xsup + ldxx;
			linalg::scale (2 * nrhs * (n + 4), 0.0, x);
		
			if (id != 0) {
				linalg::copy (nrhs, sub, xsub, lda, ldx);
				linalg::copy (nrhs, diag, xdiag, lda, ldx);
				linalg::scale (nrhs, 0.0, xsup, ldx);
			}
			if (id != np - 1) {
				linalg::scale (nrhs, 0.0, xsub + 1, ldx);
				linalg::copy (nrhs, diag + n + ntop + nbot - 1, xdiag + 1, lda, ldx);
				linalg::copy (nrhs, sup + n + ntop + nbot - 1, xsup + 1, lda, ldx);
			}
		
			for (int i = 0; i < nrhs; ++i) {
				linalg::tridiagonal_factorize (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, ipiv + i * lda, info);
			}
		
	#ifdef _MPI
			if (id != 0) {
				linalg::copy (nrhs, sub + 1, bufferl, lda, n);
				for (int i = 0; i < nrhs; ++i) {
					linalg::tridiagonal_solve (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, ipiv + i * lda, &bufferl [i * n], info);
				}
			}
		
			if (id != np - 1) {
				linalg::copy (nrhs, sup + n + ntop + nbot - 2, bufferr + n - 1, lda, n);
				for (int i = 0; i < nrhs; ++i) {
					linalg::tridiagonal_solve (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, ipiv + i * lda, &bufferr [i * n], info);
				}
			}
			if (id != 0) {
				for (int i = 0; i < nrhs; ++i) {
					xdiag [i * ldx] -= bufferl [i * n] * sup [i * lda];
				}
				if (id != np - 1) {
					for (int i = 0; i < nrhs; ++i) {
						xsup [i * ldx] -= bufferr [i * n] * sup [i * lda];
					}
				}
			}
			if (id != np - 1) {
				for (int i = 0; i < nrhs; ++i) {
					xdiag [1 + i * ldx] -= bufferr [n - 1 + i * n] * sub [n + ntop + nbot - 1 + i * lda];
				}
				if (id != 0) {
					for (int i = 0; i < nrhs; ++i) {
						xsub [1 + i * ldx] -= bufferl [n - 1 + i * n] * sub [n + ntop + nbot - 1 + i * lda];
					}
				}
			}
		
			if (id == 0) {
				MPI::COMM_WORLD.Gather (xsub, 6 * nrhs, MPI::DOUBLE, xsub, 6 * nrhs, MPI::DOUBLE, 0);
				int ldbuff = 2 * np;
				ldx = 6 * nrhs;
				std::vector <double> buffer (6 * np * nrhs);
				double *bsub = &buffer [0], *bdiag = bsub + 2 * np * nrhs, *bsup = bdiag + 2 * np * nrhs;
		
				for (int i = 0; i < nrhs; ++i) {
					for (int j = 0; j < np; ++j) {
						bsub [2 * j + i * ldbuff] = xsub [2 * i + j * ldx];
						bsub [2 * j + 1 + i * ldbuff] = xsub [2 * i + 1 + j * ldx];
					
						bdiag [2 * j + i * ldbuff] = xdiag [2 * i + j * ldx];
						bdiag [2 * j + 1 + i * ldbuff] = xdiag [2 * i + 1 + j * ldx];
					
						bsup [2 * j + i * ldbuff] = xsup [2 * i + j * ldx];
						bsup [2 * j + 1 + i * ldbuff] = xsup [2 * i + 1 + j * ldx];
					}
				}
			
				xdiag = xsub + 2 * nrhs * np;
				xsup = xdiag + 2 * nrhs * np;
				xsupsup = xsup + 2 * nrhs * np;
				linalg::copy (2 * nrhs * np, &bsub [0], xsub);
				linalg::copy (2 * nrhs * np, &bdiag [0], xdiag);
				linalg::copy (2 * nrhs * np, &bsup [0], xsup);
			
		
				for (int i = 0; i < nrhs; ++i) {
					linalg::tridiagonal_factorize (2 * np - 2, xsub + 2 + i * ldbuff, xdiag + 1 + i * ldbuff, xsup + 1 + i * ldbuff, xsupsup + 1 + i * ldbuff, xipiv + i * ldbuff, info);
				}
			
			} else {
				MPI::COMM_WORLD.Gather (xsub, 6 * nrhs, MPI::DOUBLE, xsub, 6 * nrhs, MPI::DOUBLE, 0);
			}
	#endif
		}
	
		void tridiag_solve (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, int* ipiv, double* b, double *x, int *xipiv, int *info, int nrhs, int lda, int ldb) {
			int ntop, nbot;
			std::vector <double> y (2 * np * nrhs, 0.0);
			if (id == 0) {
				ntop = 0;
			} else {
				ntop = 1;
			}
			if (id == np - 1) {
				nbot = 0;
			} else {
				nbot = 1;
			}
		
			if (lda == -1) {
				lda = n + ntop + nbot;
			}
			if (ldb == -1) {
				ldb = n + ntop + nbot;
			}
			int ldy = 2;
		
			int ldxx = 2 * nrhs;
			if (id == 0) {
				ldxx *= np;
			}
			double *bufferl = x, *bufferr = bufferl + n * nrhs, *xsub = bufferr + n * nrhs, *xdiag = xsub + ldxx, *xsup = xdiag + ldxx, *xsupsup = xsup + ldxx;
		
			linalg::add_scaled (ntop * nrhs, 1.0, b, &y [0], ldb, ldy);
			linalg::add_scaled (nbot * nrhs, 1.0, b + ntop + n, &y [1], ldb, ldy);
		
			for (int i = 0; i < nrhs; ++i) {
				linalg::tridiagonal_solve (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, ipiv + i * lda, b + ntop + i * ldb, info);
			}
		
	#ifdef _MPI
			if (id != 0) {
				for (int i = 0; i < nrhs; ++i) {
					y [i * ldy] -= sup [i * lda] * b [1 + i * ldb];
				}
			}
			if (id != np - 1) {
				for (int i = 0; i < nrhs; ++i) {
					y [1 + i * ldy] -= sub [n + ntop + nbot - 1 + i * lda] * b [n + ntop + nbot - 2 + i * ldb];
				}
			}
		
			if (id == 0) {
				MPI::COMM_WORLD.Gather (&y [0], 2 * nrhs, MPI::DOUBLE, &y [0], 2 * nrhs, MPI::DOUBLE, 0);
				int ldx = 2 * np;
				ldy = 2 * nrhs;
				int ldbuff = 2 * np;
			
				std::vector <double> buffer (2 * np * nrhs);
			
				for (int i = 0; i < nrhs; ++i) {
					for (int j = 0; j < np; ++j) {
						buffer [2 * j + i * ldbuff] = y [2 * i + j * ldy];
						buffer [2 * j + 1 + i * ldbuff] = y [2 * i + 1 + j * ldy];
					}
				}

				linalg::copy (2 * nrhs * np, &buffer [0], &y [0]);
			
				for (int i = 0; i < nrhs; ++i) {
					linalg::tridiagonal_solve (2 * np - 2, xsub + 2 + i * ldx, xdiag + 1 + i * ldx, xsup + 1 + i * ldx, xsupsup + 1 + i * ldx, xipiv + i * ldx, &y [1 + i * ldbuff], info);
				}
			
				for (int i = 0; i < nrhs; ++i) {
					for (int j = 0; j < np; ++j) {
						buffer [2 * i + j * ldy] = y [2 * j + i * ldbuff];
						buffer [2 * i + 1 + j * ldy] = y [2 * j + 1 + i * ldbuff];
					}
				}

				linalg::copy (2 * nrhs * np, &buffer [0], &y [0]);
			
				MPI::COMM_WORLD.Scatter (&y [0], 2 * nrhs, MPI::DOUBLE, &y [0], 2 * nrhs, MPI::DOUBLE, 0);
			} else {
				MPI::COMM_WORLD.Gather (&y [0], 2 * nrhs, MPI::DOUBLE, NULL, 2 * nrhs, MPI::DOUBLE, 0);
				MPI::COMM_WORLD.Scatter (NULL, 2 * nrhs, MPI::DOUBLE, &y [0], 2 * nrhs, MPI::DOUBLE, 0);
			}
			ldy = 2;
				
			if (id != 0) {
				linalg::copy (nrhs, &y [0], b, ldy, ldb);
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j < nrhs; ++j) {
						b [ntop + i + j * ldb] -= bufferl [i + j * n] * y [j * ldy]; 
					}
				}
			}
			if (id != np - 1) {
				linalg::copy (nrhs, &y [1], b + n + ntop + nbot - 1, ldy, ldb);
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j < nrhs; ++j) {
						b [ntop + i + j * ldb] -= bufferr [i + j * n] * y [1 + j * ldy];
					}
				}
			}
	#endif
		}
	} /* block */
} /* linalg */

#endif /* end of include guard: TRIDIAGONAL_CPP_53CA7D2B */
