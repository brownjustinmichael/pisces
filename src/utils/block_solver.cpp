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
	
	void p_block_tridiag_factorize (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, int* ipiv, double *x, int *xipiv, int *info, int nrhs, int lda) {
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
			tridiagonal_factorize (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, ipiv + i * lda, info);
		}
		
#ifdef _MPI
		if (id != 0) {
			copy (nrhs, sub + 1, bufferl, lda, n);
			for (int i = 0; i < nrhs; ++i) {
				tridiagonal_solve (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, ipiv + i * lda, &bufferl [i * n], info);
			}
		}
		
		if (id != np - 1) {
			copy (nrhs, sup + n + ntop + nbot - 2, bufferr + n - 1, lda, n);
			for (int i = 0; i < nrhs; ++i) {
				tridiagonal_solve (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, ipiv + i * lda, &bufferr [i * n], info);
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
				tridiagonal_factorize (2 * np - 2, xsub + 2 + i * ldbuff, xdiag + 1 + i * ldbuff, xsup + 1 + i * ldbuff, xsupsup + 1 + i * ldbuff, xipiv + i * ldbuff, info);
			}
			
		} else {
			MPI::COMM_WORLD.Gather (xsub, 6 * nrhs, MPI::DOUBLE, xsub, 6 * nrhs, MPI::DOUBLE, 0);
		}
#endif
	}
	
	void p_block_tridiag_solve (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, int* ipiv, double* b, double *x, int *xipiv, int *info, int nrhs, int lda, int ldb) {
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
		
		add_scaled (ntop * nrhs, 1.0, b, &y [0], ldb, ldy);
		add_scaled (nbot * nrhs, 1.0, b + ntop + n, &y [1], ldb, ldy);
		
		for (int i = 0; i < nrhs; ++i) {
			tridiagonal_solve (n, sub + 1 + ntop + i * lda, diag + ntop + i * lda, sup + ntop + i * lda, supsup + ntop + i * lda, ipiv + i * lda, b + ntop + i * ldb, info);
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

			copy (2 * nrhs * np, &buffer [0], &y [0]);
			
			for (int i = 0; i < nrhs; ++i) {
				tridiagonal_solve (2 * np - 2, xsub + 2 + i * ldx, xdiag + 1 + i * ldx, xsup + 1 + i * ldx, xsupsup + 1 + i * ldx, xipiv + i * ldx, &y [1 + i * ldbuff], info);
			}
			
			for (int i = 0; i < nrhs; ++i) {
				for (int j = 0; j < np; ++j) {
					buffer [2 * i + j * ldy] = y [2 * j + i * ldbuff];
					buffer [2 * i + 1 + j * ldy] = y [2 * j + 1 + i * ldbuff];
				}
			}

			copy (2 * nrhs * np, &buffer [0], &y [0]);
			
			MPI::COMM_WORLD.Scatter (&y [0], 2 * nrhs, MPI::DOUBLE, &y [0], 2 * nrhs, MPI::DOUBLE, 0);
		} else {
			MPI::COMM_WORLD.Gather (&y [0], 2 * nrhs, MPI::DOUBLE, NULL, 2 * nrhs, MPI::DOUBLE, 0);
			MPI::COMM_WORLD.Scatter (NULL, 2 * nrhs, MPI::DOUBLE, &y [0], 2 * nrhs, MPI::DOUBLE, 0);
		}
		ldy = 2;
				
		if (id != 0) {
			copy (nrhs, &y [0], b, ldy, ldb);
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < nrhs; ++j) {
					b [ntop + i + j * ldb] -= bufferl [i + j * n] * y [j * ldy]; 
				}
			}
		}
		if (id != np - 1) {
			copy (nrhs, &y [1], b + n + ntop + nbot - 1, ldy, ldb);
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < nrhs; ++j) {
					b [ntop + i + j * ldb] -= bufferr [i + j * n] * y [1 + j * ldy];
				}
			}
		}
#endif
	}
	
	void p_block_banded_factorize (int id, int np, int n, int kl, int ku, double* matrix, int* ipiv, double *x, int *xipiv, double *bufferl, double *bufferr, int *info, int nrhs, int lda) {
		int ntop, nbot;
		std::stringstream debug;
		int ldx = (kl + ku);
		if (id == 0) {
			ntop = 0;
		} else {
			ntop = ku;
		}
		if (id == np - 1) {
			nbot = 0;
		} else {
			nbot = kl;
		}

		if (lda == -1) {
			lda = 2 * kl + ku + 1;
		}
		int ldaa = n + ntop + nbot;

		int ldxx = 2 * ku + 2 * kl;
		utils::scale (nrhs * ldx, 0.0, x);
		utils::scale (nrhs * kl, 0.0, bufferl);
		utils::scale (nrhs * ku, 0.0, bufferl);

		for (int i = 0; i < nrhs; ++i) {
			matrix_banded_factorize (n, n, kl, ku, matrix + i * lda * ldaa + ntop * lda, ipiv + i * n, info, lda);
		}

#ifdef _MPI
		if (id != 0) {
			for (int i = 0; i < nrhs; ++i) {
				for (int j = 0; j < kl + ku; ++j) {
					for (int k = 0; k < ku; ++k) {
						x [i * ldxx * ldx + j * ldx + k] = matrix [i * ldaa * lda + j * lda + kl + ku + 1 - j];
					}
				}
				for (int j = 0; j < kl; ++j) {
					for (int k = 0; k < j; ++k) {
						bufferl [i * n * kl + j * n + k] = matrix [i * lda * ldaa + (j + ku) * lda + k + 2 * kl + ku - j];
					}
				}
				
				for (int k = 0; k < n; ++k) {
					for (int j = 0; j < kl; ++j) {
						debug << bufferl [i * nrhs * kl + j * n + k] << " "; 
					}
					DEBUG (debug.str ());
					debug.str ("");
				}
				
				matrix_banded_solve (n, kl, ku, matrix + i * lda * ldaa + ntop * lda, ipiv + i * n, &bufferl [i * n * kl], info);
				for (int j = 0; j < ku; ++j) {
					for (int k = 0; k < ku; ++k) {
						for (int l = 0; l < j + 1; ++l) {
							x [i * ldxx * ldx + (j + ku) * ldx + k] -= bufferl [i * n * kl + j * n + l] * matrix [i * lda * ldaa + (kl + ku + l) * lda + kl + ku - l];
						}
					}
					if (id != np - 1) {
						for (int k = 0; k < kl; ++k) {
							for (int l = 0; l < j + 1; ++l) {
								x [i * ldxx * ldx + (j + kl) * ldx + k] -= bufferr [i * n * ku + j * n + l] * matrix [i * lda * ldaa + (j + ldaa - ku - 2 * kl + l) * lda + kl + ku - l];
							}
						}
					}
				}
			}
		}

		if (id != np - 1) {
			for (int i = 0; i < nrhs; ++i) {
				for (int j = 0; j < ku + kl; ++j) {
					for (int k = 0; k < kl; ++k) {
						x [i * ldxx * ldx + (j + ku + kl) * ldx + ku + 1 + k] = matrix [i * ldaa * lda + (j + ldaa - ku - kl) * lda + kl + ku + k - j];
					}
				}
				for (int j = 0; j < ku; ++j) {
					for (int k = -1; k >= -ku; --k) {
						bufferr [i * n * ku + j * n + k + n] = matrix [i * lda * ldaa + (ldaa - kl - ku + j) * lda + k + kl + ku - j];
					}
				}
				matrix_banded_solve (n, kl, ku, matrix + i * lda * ldaa + ntop * lda, ipiv + i * n, &bufferr [i * n * ku], info);
				for (int j = 0; j < kl; ++j) {
					for (int k = 0; k < kl; ++k) {
						for (int l = 0; l < j + 1; ++l) {
							x [i * ldxx * ldx + (j + ku + kl) * ldx + k] -= bufferl [i * n * kl + j * n + l] * matrix [i * ldaa * lda + (l + ldaa - ku - 2 * kl) * lda + kl + ku + k - l];
						}
					}
					if (id != 0) {
						for (int k = 0; k < kl; ++k) {
							for (int l = 0; l < j + 1; ++l) {
								x [i * ldxx * ldx + (j + ku + kl) * ldx + k] -= bufferr [i * n * ku + j * n + l] * matrix [i * lda * ldaa + (kl + ku + l) * lda + kl + ku - l];
							}
						}
					}
				}
			}
		}

		DEBUG ("GATHER");

		if (id == 0) {
			ldx = np * (2 * ku + 2 * kl);
			std::vector <double> buffer (np * (2 * ku + 2 * kl) * (ku + kl));
			for (int q = 0; q < nrhs; ++q) {
				MPI::COMM_WORLD.Gather (x + q * ldx * ldx, (2 * ku + 2 * kl) * (ku + kl), MPI::DOUBLE, &buffer [0], (2 * ku + 2 * kl) * (ku + kl), MPI::DOUBLE, 0);
				
				for (int i = 0; i < (np - 1) * (kl + ku); ++i) {
					scale ((np - 1) * (kl + ku), 0.0, x + q * ldx * ldx + i * ldx);
				}
				int bcur = 0, xcur = 0;
				for (int i = 0; i < np - 1; ++i) {
					for (int j = 0; j < 2 * ku + 2 * kl; ++j) {
						for (int k = kl; k < 2 * kl + ku; ++k) {
							x [q * ldx * ldx + xcur + k * ldx + j] += buffer [bcur + (k - kl) * (kl + ku) + j];
						}
					}
					bcur += (2 * ku + 2 * kl) * (ku + kl);
					xcur += (ku + kl) * (ldx + 1);
				}
				for (int j = 0; j < 2 * ku + 2 * kl; ++j) {
					for (int k = kl; k < 2 * kl + ku; ++k) {
						x [q * ldx * ldx + xcur + k * ldx + j] += buffer [bcur + (k - kl) * (kl + ku) + j];
					}
				}
				
				matrix_factorize ((ku + kl) * (np - 1), (ku + kl) * (np - 1), x + q * ldx * ldx, xipiv + q * ldx, info, ldx);
			}

		} else {
			MPI::COMM_WORLD.Gather (x, (2 * ku + 2 * kl) * (ku + kl), mpi_type (&typeid (double)), NULL, (2 * ku + 2 * kl) * (ku + kl), mpi_type (&typeid (double)), 0);
		}
#endif
	}

	void p_block_banded_solve (int id, int np, int n, int kl, int ku, double* matrix, int* ipiv, double* b, double *x, int *xipiv, double *bufferl, double *bufferr, int *info, int nrhs, int lda, int ldb) {
		DEBUG ("BEGIN");
		std::stringstream debug;
		std::vector <double> y (2 * (ku + kl) * np * nrhs, 0.0);
		DEBUG ("y len " << 2 * (ku + kl) * np * nrhs);
		int ntop, nbot;
		int ldx = (kl + ku);
		if (id == 0) {
			ntop = 0;
		} else {
			ntop = ku;
		}
		if (id == np - 1) {
			nbot = 0;
		} else {
			nbot = kl;
		}
		
		int ldy = 0;
		if (id == 0) {
			ldy = 2 * (kl + ku) * np;
		} else {
			ldy = 2 * (kl + ku);
		}

		if (lda == -1) {
			lda = 2 * kl + ku + 1;
		}
		if (ldb == -1) {
			ldb = n + ntop + nbot;
		}
		int ldaa = n + ntop + nbot;

		int ldxx = 2 * ku + 2 * kl;

		add_scaled (ntop * nrhs, 1.0, b, &y [kl], ldb, ldy);
		add_scaled (nbot * nrhs, 1.0, b + ntop + n, &y [2 * kl - ku], ldb, ldy);

		for (int i = 0; i < nrhs; ++i) {
			matrix_banded_solve (n, kl, ku, matrix + i * lda * ldaa + ntop * lda, ipiv + i * n, b + ntop + i * ldb, info);
		}
		

// #ifdef _MPI
// 		if (id != 0) {
// 			for (int j = 0; j < nrhs; ++j) {
// 				for (int k = 0; k < ku; ++k) {
// 					for (int l = 0; l < j + 1; ++l) {
// 						y [j * ldy + k + kl] -= b [j * ldb + l + ku] * matrix [j * lda * ldaa + (kl + ku + l) * lda + kl + ku - l];
// 					}
// 				}
// 			}
// 		}
// 		if (id != np - 1) {
// 			for (int j = 0; j < nrhs; ++j) {
// 				for (int k = 0; k < kl; ++k) {
// 					for (int l = 0; l < j + 1; ++l) {
// 						y [j * ldx + k + kl + ku] -= b [j * ldb + l + ku] * matrix [j * ldaa * lda + (l + ldaa - ku - 2 * kl) * lda + kl + ku + k - l];
// 					}
// 				}
// 			}
// 		}
// 		if (id == 0) {
// 			int ycur = 0, bcur = 0;
// 			// DEBUG ("GATHERING" << &y [0]);
// 			std::vector <double> buffer (2 * ldy);
// 			// DEBUG ("BUFFER " << &buffer [0]);
// 			for (int j = 0; j < nrhs; ++j) {
// 				MPI::COMM_WORLD.Gather (&y [j * ldy], 2 * (kl + ku), mpi_type (&typeid (double)), &buffer [0], 2 * (kl + ku), mpi_type (&typeid (double)), 0);
// 				utils::scale (ldy, 0.0, &y [0]);
// 				for (int i = 0; i < np; ++i) {
// 						for (int k = 0; k < 2 * (ku + kl); ++k) {
// 							// DEBUG ("RUN " << j * ldy + ycur + k << " " << bcur + k << " " << buffer [bcur + k]);
// 							y [j * ldy + ycur + k] += buffer [bcur + k];
// 						}
// 					ycur += ku + kl;
// 					bcur += (ku + kl);
// 				}
//
// 				// for (int i = 0; i < ldy; ++i) {
// // 					DEBUG ("Y = " << y [j * ldy + i] << " " << j * ldy + i);
// // 				}
//
// 				// DEBUG ("Solving");
// 				matrix_solve ((ku + kl) * (np - 1), x + j * ldx * ldx, xipiv + j * ldx, &y [j * ldy + ku], info, 1, ldx);
// 				// DEBUG ("Solved");
//
// 				ycur = 0;
// 				bcur = 0;
// 				for (int i = 0; i < np; ++i) {
// 					for (int k = 0; k < (ku + kl); ++k) {
// 						buffer [bcur + k] = y [j * ldy + ycur + k];
// 					}
// 					ycur += ku + kl;
// 					bcur += (ku + kl);
// 				}
//
// 				MPI::COMM_WORLD.Scatter (&buffer [0], 2 * (kl + ku), mpi_type (&typeid (double)), &y [j * ldy], 2 * (kl + ku), mpi_type (&typeid (double)), 0);
// 			}
// 		} else {
// 			for (int j = 0; j < nrhs; ++j) {
// 				MPI::COMM_WORLD.Gather (&y [j * ldy], 2 * (kl + ku), mpi_type (&typeid (double)), NULL, 2 * (kl + ku), mpi_type (&typeid (double)), 0);
// 				MPI::COMM_WORLD.Scatter (NULL, 2 * (kl + ku), mpi_type (&typeid (double)), &y [j * ldy], 2 * (kl + ku), mpi_type (&typeid (double)), 0);
// 			}
// 		}
// 		DEBUG ("OUT");
// 		matrix_copy (ntop, nrhs, &y [0], b, ntop + nbot, ldb);
// 		matrix_copy (nbot, nrhs, &y [ntop], b + ntop + n, ntop + nbot, ldb);
//
// 		matrix_matrix_multiply (n, nrhs, ntop, -1.0, &bufferl [0], b, 1.0, b + ntop, lda, ldb, ldb);
// 		matrix_matrix_multiply (n, nrhs, nbot, -1.0, &bufferr [0], b + ntop + n, 1.0, b + ntop, lda, ldb, ldb);
// 		DEBUG ("DONE");
// #endif
	}
} /* utils */