/*!**********************************************************************
 * \file block_solver.cpp
 * /Users/justinbrown/Dropbox/pisces/src/utils
 * 
 * Created by Justin Brown on 2013-11-17.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <vector>
#ifdef _MPI
#include <mpi.h>
#endif
#include "../utils/utils.hpp"
#include "../utils/solver_utils.hpp"
#include "block_solver.hpp"

namespace utils
{
// 	void p_block_matrix_factorize (int id, int np, int n, int ntop, int nbot, double* a, int* ipiv, double *x, int *xipiv, int* ns, int *info, int lda, int ldx) {
// 		int ntot = 0, n2tot = 0, ldb = ntop + nbot;
// 		double* am = a + ntop * (lda + 1);
// 		double* at = a + ntop * lda;
// 		double* ab = a + ntop * lda + ntop + n;
// 		double* ar = a + (ntop + n) * lda + ntop;
// 		double* al = a + ntop;
// 		if (lda == -1) {
// 			lda = n + ntop + nbot;
// 		}
// 		
// 		std::vector <int> ns2;
// 		if (id == 0) {
// 			ns2.resize (np);
// 			for (int i = 0; i < np - 1; ++i) {
// 				ntot += ns [i];
// 				ns2 [i] = (ns [i] + ns [i + 1]) * (ns [i] + ns [i + 1]);
// 				n2tot += ns2 [i];
// 			}
// 			ntot += ns [np - 1];
// 			ns2 [np - 1] = ns [np - 1] * ns [np - 1];
// 			n2tot += ns2 [np - 1];
// 		} else {
// 			n2tot = ntop + nbot;
// 		}
// 						
// 		if (ldx == -1) {
// 			if (id == 0) {
// 				ldx = ntot;
// 			} else {
// 				ldx = ntop + nbot;
// 			}
// 		}
// 		if (ldx < ntop + nbot) {
// 			throw 0;
// 		}
// 				
// 		matrix_copy (ntop, ntop, a, x, lda, ldb);
// 		matrix_copy (nbot, ntop, a + ntop + n, x + ntop, lda, ldb);
// 		matrix_copy (ntop, nbot, a + (ntop + n) * lda, x + ntop * ldb, lda, ldb);
// 		matrix_copy (nbot, nbot, a + (ntop + n) * (lda + 1), x + ntop * (ldb + 1), lda, ldb);
// 		
// 		matrix_factorize (n, n, am, ipiv, info, lda);
// 		
// 		matrix_solve (n, am, ipiv, al, info, ntop, lda, lda);
// 		matrix_solve (n, am, ipiv, ar, info, nbot, lda, lda);
// 		
// #ifdef _MPI
// 		matrix_matrix_multiply (ntop, ntop, n, -1.0, at, al, 1.0, x, lda, lda, ldb);
// 		matrix_matrix_multiply (nbot, ntop, n, -1.0, ab, al, 1.0, x + ntop, lda, lda, ldb);
// 		matrix_matrix_multiply (ntop, nbot, n, -1.0, at, ar, 1.0, x + ntop * ldx, lda, lda, ldb);
// 		matrix_matrix_multiply (nbot, nbot, n, -1.0, ab, ar, 1.0, x + ntop * (ldx + 1), lda, lda, ldb);
// 
// 		if (id == 0) {
// 			std::vector <double> buffer (n2tot);
// 			std::vector <int> displs;
// 			if (id == 0) {
// 				displs.resize (np);
// 				for (int i = 1; i < np; ++i) {
// 					displs [i] = displs [i - 1] + ns2 [i - 1];
// 				}
// 			}
// 			MPI::COMM_WORLD.Gatherv (x, (ntop + nbot) * (ntop + nbot), MPI::DOUBLE, &buffer [0], &ns2 [0], &displs [0], MPI::DOUBLE, 0);
// 
// 			for (int i = 0; i < ntot; ++i) {
// 				scale (ntot, 0.0, x + i * ldx);
// 			}
// 			int bcur = 0, xcur = 0;
// 			for (int i = 0; i < np - 1; ++i) {
// 				for (int j = 0; j < ns [i] + ns [i + 1]; ++j) {
// 					for (int k = 0; k < ns [i] + ns [i + 1]; ++k) {
// 						x [xcur + k * ldx + j] += buffer [bcur + k * (ns [i] + ns [i + 1]) + j];
// 					}
// 				}
// 				bcur += (ns [i] + ns [i + 1]) * (ns [i] + ns [i + 1]);
// 				xcur += ns [i] * (ldx + 1);
// 			}
// 			for (int j = 0; j < ns [np - 1]; ++j) {
// 				for (int k = 0; k < ns [np - 1]; ++k) {
// 					x [xcur + k * ldx + j] += buffer [bcur + k * ns [np - 1] + j];
// 				}
// 			}
// 			
// 			matrix_factorize (ntot, ntot, x, xipiv, info, ldx);
// 		} else {
// 			MPI::COMM_WORLD.Gatherv (x, (ntop + nbot) * (ntop + nbot), MPI::DOUBLE, NULL, NULL, NULL, MPI::DOUBLE, 0);
// 		}
// #endif
// 		/*
// 			TODO Check that these work w/o MPI
// 		*/
// 	}
	
// 	void p_block_matrix_solve (int id, int np, int n, int ntop, int nbot, double* a, int* ipiv, double* b, double *x, int *xipiv, int *ns, int *info, int nrhs, int lda, int ldx, int ldb) {
// 		std::vector <int> nsp (np);
// 		double* am = a + ntop * (lda + 1);
// 		double* at = a + ntop * lda;
// 		double* ab = a + ntop * lda + ntop + n;
// 		double* ar = a + (ntop + n) * lda + ntop;
// 		double* al = a + ntop;
// 		
// 		int ldy = 0;
// 		if (id == 0) {
// 			for (int i = 0; i < np - 1; ++i) {
// 				ldy += ns [i];
// 				nsp [i] = (ns [i] + ns [i + 1]) * nrhs;
// 			}
// 			ldy += ns [np - 1];
// 			nsp [np - 1] = (ns [np - 1]) * nrhs;
// 		} else {
// 			ldy = ntop + nbot;
// 		}
// 		
// 		std::vector <double> y (ldy * nrhs);
// 		
// 		if (lda == -1) {
// 			lda = n + ntop + nbot;
// 		}
// 		if (ldb == -1) {
// 			ldb = n + ntop + nbot;
// 		}
// 		if (ldx == -1) {
// 			ldx = ldy;
// 		}
// 
// 		matrix_add_scaled (ntop, nrhs, 1.0, b, &y [0], ldb, ntop + nbot);
// 		matrix_add_scaled (nbot, nrhs, 1.0, b + ntop + n, &y [ntop], ldb, ntop + nbot);
// 		
// 		matrix_solve (n, am, ipiv, b + ntop, info, nrhs, lda, ldb);
// 				
// #ifdef _MPI
// 		matrix_matrix_multiply (ntop, nrhs, n, -1.0, at, b + ntop, 1.0, &y [0], lda, ldb, ntop + nbot);
// 		matrix_matrix_multiply (nbot, nrhs, n, -1.0, ab, b + ntop, 1.0, &y [ntop], lda, ldb, ntop + nbot);
// 			
// 		if (id == 0) {
// 			int ycur = 0, bcur = 0;
// 			std::vector <double> buffer (nrhs * 2 * ldy);
// 			std::vector <int> displs;
// 			if (id == 0) {
// 				displs.resize (np);
// 				for (int i = 1; i < np; ++i) {
// 					displs [i] = displs [i - 1] + nsp [i - 1];
// 				}
// 			}
// 			MPI::COMM_WORLD.Gatherv (&y [0], (ntop + nbot) * nrhs, MPI::DOUBLE, &buffer [0], &nsp [0], &displs [0], MPI::DOUBLE, 0);
// 			utils::scale (nrhs * ldy, 0.0, &y [0]);
// 			for (int i = 0; i < np; ++i) {
// 				for (int j = 0; j < nrhs; ++j) {
// 					for (int k = 0; k < nsp [i] / nrhs; ++k) {
// 						y [j * ldy + ycur + k] += buffer [bcur + j * nsp [i] / nrhs + k];
// 					}
// 				}
// 				ycur += ns [i];
// 				bcur += nsp [i];
// 			}
// 			
// 			matrix_solve (ldy, x, xipiv, &y [0], info, nrhs, ldx, ldy);
// 			
// 			ycur = 0;
// 			bcur = 0;
// 			for (int i = 0; i < np; ++i) {
// 				for (int j = 0; j < nrhs; ++j) {
// 					for (int k = 0; k < nsp [i] / nrhs; ++k) {
// 						buffer [bcur + j * nsp [i] / nrhs + k] = y [j * ldy + ycur + k];
// 					}
// 				}
// 				ycur += ns [i];
// 				bcur += nsp [i];
// 			}
// 			
// 			MPI::COMM_WORLD.Scatterv (&buffer [0], &nsp [0], &displs [0], MPI::DOUBLE, &y [0], (ntop + nbot) * nrhs, MPI::DOUBLE, 0);
// 		} else {
// 			MPI::COMM_WORLD.Gatherv (&y [0], (ntop + nbot) * nrhs, MPI::DOUBLE, NULL, NULL, NULL, MPI::DOUBLE, 0);
// 			MPI::COMM_WORLD.Scatterv (NULL, NULL, NULL, MPI::DOUBLE, &y [0], (ntop + nbot) * nrhs, MPI::DOUBLE, 0);
// 		}
// 		
// 		matrix_copy (ntop, nrhs, &y [0], b, ntop + nbot, ldb);
// 		matrix_copy (nbot, nrhs, &y [ntop], b + ntop + n, ntop + nbot, ldb);
// 		
// 		matrix_matrix_multiply (n, nrhs, ntop, -1.0, al, b, 1.0, b + ntop, lda, ldb, ldb);
// 		matrix_matrix_multiply (n, nrhs, nbot, -1.0, ar, b + ntop + n, 1.0, b + ntop, lda, ldb, ldb);
// #endif
// 	}
	
	void p_block_tridiag_factorize (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, int* ipiv, double *x, int *xipiv, int *info, int nrhs, int lda, int inrhs) {
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
		utils::scale (2 * nrhs * (n + 4), 0.0, x);
		
		if (id != 0) {
			copy (nrhs, sub, xsub, lda, ldx);
			copy (nrhs, diag, xdiag, lda, ldx);
			scale (nrhs, 0.0, xsup, ldx);
		}
		if (id != np - 1) {
			scale (nrhs, 0.0, xsub + 1, ldx);
			copy (nrhs, diag + n + ntop + nbot - 1, xdiag + 1, lda, ldx);
			copy (nrhs, sup + n + ntop + nbot - 1, xsup + 1, lda, ldx);
		}
		
		for (int i = 0; i < nrhs; ++i) {
			// tridiagonal_factorize (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, ipiv + i * lda, info);
		}
		
#ifdef _MPI
		if (id != 0) {
			copy (nrhs, sub + 1, bufferl, lda, n);
			for (int i = 0; i < nrhs; ++i) {
				// tridiagonal_solve (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, ipiv + i * lda, &bufferl [i * n], info);
				tridiagonal_direct_solve (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, &bufferl [i * n]);
			}
		}
		
		if (id != np - 1) {
			copy (nrhs, sup + n + ntop + nbot - 2, bufferr + n - 1, lda, n);
			for (int i = 0; i < nrhs; ++i) {
				// tridiagonal_solve (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, ipiv + i * lda, &bufferr [i * n], info);
				tridiagonal_direct_solve (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, &bufferr [i * n]);
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
			copy (2 * nrhs * np, &bsub [0], xsub);
			copy (2 * nrhs * np, &bdiag [0], xdiag);
			copy (2 * nrhs * np, &bsup [0], xsup);
			
		
			for (int i = 0; i < nrhs; ++i) {
				// tridiagonal_factorize (2 * np - 2, xsub + 2 + i * ldbuff, xdiag + 1 + i * ldbuff, xsup + 1 + i * ldbuff, xsupsup + 1 + i * ldbuff, xipiv + i * ldbuff, info);
			}
			
		} else {
			MPI::COMM_WORLD.Gather (xsub, 6 * nrhs, MPI::DOUBLE, xsub, 6 * nrhs, MPI::DOUBLE, 0);
		}
#endif
	}
	
	void p_block_tridiag_solve (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, int* ipiv, double* b, double *x, int *xipiv, int *info, int nrhs, int lda, int ldb, int inrhs) {
		int ntop, nbot;
		std::vector <double> y (2 * np * nrhs * inrhs, 0.0);
		
		for (int i = 0; i < inrhs; ++i) {
			for (int j = 0; j < nrhs; ++j) {
				for (int k = 0; k < n; ++k) {
					if (b [k + j * ldb + k * nrhs * ldb] != b [k + j * ldb + k * nrhs * ldb]) {
						DEBUG ("Found nan at beginning");
						throw 0;
					}
				}
			}
		}
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
		
		for (int j = 0; j < inrhs; ++j) {
			add_scaled (ntop * nrhs, 1.0, b + j * nrhs * ldb, &y [j * nrhs * ldy], ldb, ldy);
			add_scaled (nbot * nrhs, 1.0, b + ntop + n + j * nrhs * ldb, &y [1 + j * nrhs * ldy], ldb, ldy);
		}
		
		for (int i = 0; i < inrhs; ++i) {
			for (int j = 0; j < nrhs; ++j) {
				for (int k = 0; k < n; ++k) {
					if (b [k + j * ldb + k * nrhs * ldb] != b [k + j * ldb + k * nrhs * ldb]) {
						DEBUG ("Found nan after scale");
						throw 0;
					}
				}
			}
		}
		
		for (int i = 0; i < nrhs; ++i) {
			// tridiagonal_solve (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, ipiv + i * lda, b + ntop + i * ldb, info, inrhs, ldb * nrhs);
			tridiagonal_direct_solve (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, b + ntop + i * ldb, inrhs, ldb * nrhs);
		}
		
		for (int i = 0; i < inrhs; ++i) {
			for (int j = 0; j < nrhs; ++j) {
				for (int k = 0; k < n; ++k) {
					if (b [k + j * ldb + k * nrhs * ldb] != b [k + j * ldb + k * nrhs * ldb]) {
						DEBUG ("Found nan after first solve");
						throw 0;
					}
				}
			}
		}
		
#ifdef _MPI
		if (id != 0) {
			for (int j = 0; j < inrhs; ++j) {
				for (int i = 0; i < nrhs; ++i) {
					y [i * ldy + j * nrhs * ldy] -= sup [i * lda] * b [1 + i * ldb + j * nrhs * ldb];
				}
			}
		}
		if (id != np - 1) {
			for (int j = 0; j < inrhs; ++j) {
				for (int i = 0; i < nrhs; ++i) {
					y [1 + i * ldy + j * nrhs * ldy] -= sub [n + ntop + nbot - 1 + i * lda] * b [n + ntop + nbot - 2 + i * ldb + j * nrhs * ldb];
				}
			}
		}
		
		for (int i = 0; i < inrhs; ++i) {
			for (int j = 0; j < nrhs; ++j) {
				for (int k = 0; k < n; ++k) {
					if (b [k + j * ldb + k * nrhs * ldb] != b [k + j * ldb + k * nrhs * ldb]) {
						DEBUG ("Found nan before gather");
						throw 0;
					}
				}
			}
		}
		
		if (id == 0) {
			MPI::COMM_WORLD.Gather (&y [0], 2 * nrhs * inrhs, MPI::DOUBLE, &y [0], 2 * nrhs * inrhs, MPI::DOUBLE, 0);
			int ldx = 2 * np;
			ldy = 2 * nrhs;
			int ldbuff = 2 * np;
			
			std::vector <double> buffer (2 * np * nrhs * inrhs);
			
			for (int k = 0; k < inrhs; ++k) {
				for (int i = 0; i < nrhs; ++i) {
					for (int j = 0; j < np; ++j) {
						buffer [2 * j + i * ldbuff + k * nrhs * ldbuff] = y [2 * i + j * inrhs * ldy + k * ldy];
						buffer [2 * j + 1 + i * ldbuff + k * nrhs * ldbuff] = y [2 * i + 1 + j * ldy * inrhs + k * ldy];
					}
				}
			}

			copy (2 * nrhs * np * inrhs, &buffer [0], &y [0]);

			for (int i = 0; i < nrhs; ++i) {
				// tridiagonal_solve (2 * np - 2, xsub + 2 + i * ldx, xdiag + 1 + i * ldx, xsup + 1 + i * ldx, xsupsup + 1 + i * ldx, xipiv + i * ldx, &y [1 + i * ldbuff], info, inrhs, nrhs * ldbuff);
				tridiagonal_direct_solve (2 * np - 2, xsub + 2 + i * ldx, xdiag + 1 + i * ldx, xsup + 1 + i * ldx, xsupsup + 1 + i * ldx, &y [1 + i * ldbuff], inrhs, nrhs * ldbuff);
			}
			
			for (int k = 0; k < inrhs; ++k) {
				for (int i = 0; i < nrhs; ++i) {
					for (int j = 0; j < np; ++j) {
						buffer [2 * i + j * ldy * inrhs + k * ldy] = y [2 * j + i * ldbuff + k * nrhs * ldbuff];
						buffer [2 * i + 1 + j * inrhs * ldy + k * ldy] = y [2 * j + 1 + i * ldbuff + k * nrhs * ldbuff];
					}
				}
			}


			copy (2 * nrhs * np * inrhs, &buffer [0], &y [0]);
			
			MPI::COMM_WORLD.Scatter (&y [0], 2 * nrhs * inrhs, MPI::DOUBLE, &y [0], 2 * nrhs * inrhs, MPI::DOUBLE, 0);
		} else {
			MPI::COMM_WORLD.Gather (&y [0], 2 * nrhs * inrhs, MPI::DOUBLE, NULL, 2 * nrhs * inrhs, MPI::DOUBLE, 0);
			MPI::COMM_WORLD.Scatter (NULL, 2 * nrhs * inrhs, MPI::DOUBLE, &y [0], 2 * nrhs * inrhs, MPI::DOUBLE, 0);
		}
		ldy = 2;
		
		for (int i = 0; i < inrhs; ++i) {
			for (int j = 0; j < nrhs; ++j) {
				for (int k = 0; k < n; ++k) {
					if (b [k + j * ldb + k * nrhs * ldb] != b [k + j * ldb + k * nrhs * ldb]) {
						DEBUG ("Found nan after scatter");
						throw 0;
					}
				}
			}
		}
		
		if (id != 0) {
			copy (nrhs * inrhs, &y [0], b, ldy, ldb);
			for (int k = 0; k < inrhs; ++k) {
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j < nrhs; ++j) {
						b [ntop + i + j * ldb + k * nrhs * ldb] -= bufferl [i + j * n] * y [j * ldy + k * nrhs * ldy]; 
					}
				}
			}
		}
		if (id != np - 1) {
			copy (nrhs * inrhs, &y [1], b + n + ntop + nbot - 1, ldy, ldb);
			for (int k = 0; k < inrhs; ++k) {
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j < nrhs; ++j) {
						b [ntop + i + j * ldb + k * nrhs * ldb] -= bufferr [i + j * n] * y [1 + j * ldy + k * nrhs * ldy];
					}
				}
			}
		}
		
		for (int i = 0; i < inrhs; ++i) {
			for (int j = 0; j < nrhs; ++j) {
				for (int k = 0; k < n; ++k) {
					if (b [k + j * ldb + k * nrhs * ldb] != b [k + j * ldb + k * nrhs * ldb]) {
						DEBUG ("Found nan at end");
						throw 0;
					}
				}
			}
		}
#endif
	}
	
	void p_block_direct_tridiag_solve (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, double* b, int nrhs, int ldb) {
		std::vector <double> y (nrhs + 1, 0.0);
		
		int ntop = 0, nbot = 0;
		if (id != 0) {
			ntop += 1;
		}
		if (id != np - 1) {
			nbot += 1;
		}
		
		if (ldb == -1) {
			ldb = n + ntop + nbot;
		}
		
		copy (n + ntop + nbot, sup, supsup);
		if (id > 0) {
			MPI::COMM_WORLD.Recv (&y [0], nrhs + 1, MPI::DOUBLE, id - 1, 0);
			supsup [0] /= diag [0] - sub [0] * y [nrhs];
			for (int i = 0; i < nrhs; ++i) {
				b [i * ldb] = (b [i * ldb] - sub [0] * y [i]) / (diag [0] - sub [0] * y [nrhs]);
			}
		} else {
			if (diag [0] == 0.0) {
				FATAL ("Algorithm can't handle 0 in first diagonal element.");
				throw 0;
			}
			supsup [0] /= diag [0];
			for (int i = 0; i < nrhs; ++i) {
				b [i * ldb] /= diag [0];
			}
		}
		
		for (int j = 1; j < n + ntop + nbot; ++j) {
			supsup [j] /= diag [j] - sub [j] * supsup [j - 1];
			for (int i = 0; i < nrhs; ++i) {
				b [i * ldb + j] = (b [i * ldb + j] - sub [j] * b [i * ldb + j - 1]) / (diag [j] - sub [j] * supsup [j - 1]);
				if (b [i * ldb + j] != b [i * ldb + j]) {
					DEBUG ("Nan at " << i << " " << j << " " << sub [j] << " " << b [i * ldb + j - 1] << " " << diag [j] << " " << supsup [j - 1]);
					throw 0;
				}
			}
		}
		
		if (id < np - 1) {
			for (int i = 0; i < nrhs; ++i) {
				y [i] = b [i * ldb + n + ntop + nbot - 1];
			}
			y [nrhs] = supsup [n + ntop + nbot - 1];
			
			MPI::COMM_WORLD.Send (&y [0], nrhs + 1, MPI::DOUBLE, id + 1, 0);
			
			MPI::COMM_WORLD.Recv (&y [0], nrhs, MPI::DOUBLE, id + 1, 1);
			
			for (int i = 0; i < nrhs; ++i) {
				b [i * ldb + n + ntop + nbot - 1] -= supsup [n + ntop + nbot - 1] * y [i];
			}
		}
		
		for (int j = n + ntop + nbot - 2; j >= 0; --j) {
			for (int i = 0; i < nrhs; ++i) {
				b [i * ldb + j] -= supsup [j] * b [i * ldb + j + 1];
			}
		}
		
		if (id > 0) {
			for (int i = 0; i < nrhs; ++i) {
				y [i] = b [i * ldb];
			}
			
			MPI::COMM_WORLD.Send (&y [0], nrhs, MPI::DOUBLE, id - 1, 1);
		}
	}
} /* utils */