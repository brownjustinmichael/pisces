/*!**********************************************************************
 * \file block_solver.hpp
 * /Users/justinbrown/Dropbox/pisces/src/utils
 * 
 * Created by Justin Brown on 2013-11-17.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BLOCK_SOLVER_HPP_8C3ZNSDI
#define BLOCK_SOLVER_HPP_8C3ZNSDI

#include "messenger/messenger.hpp"
#include "linalg/utils.hpp"
#include "linalg/linalg.hpp"
#include <vector>

namespace utils
{
	/*!**********************************************************************
	 * \brief Factorize a tridiagonal block matrix
	 ************************************************************************/
	template <class datatype>
	void p_block_matrix_factorize (int id, int np, int n, int ntop, int nbot, datatype* a, int* ipiv, datatype *x, int *xipiv, int* ns, int *info, int lda = -1, int ldx = -1) {
		int ntot = 0, n2tot = 0, ldb = ntop + nbot;
		datatype* am = a + ntop * (lda + 1);
		datatype* at = a + ntop * lda;
		datatype* ab = a + ntop * lda + ntop + n;
		datatype* ar = a + (ntop + n) * lda + ntop;
		datatype* al = a + ntop;
		if (lda == -1) {
			lda = n + ntop + nbot;
		}
		
		std::vector <int> ns2;
		if (id == 0) {
			ns2.resize (np);
			for (int i = 0; i < np - 1; ++i) {
				ntot += ns [i];
				ns2 [i] = (ns [i] + ns [i + 1]) * (ns [i] + ns [i + 1]);
				n2tot += ns2 [i];
			}
			ntot += ns [np - 1];
			ns2 [np - 1] = ns [np - 1] * ns [np - 1];
			n2tot += ns2 [np - 1];
		} else {
			n2tot = ntop + nbot;
		}
						
		if (ldx == -1) {
			if (id == 0) {
				ldx = ntot;
			} else {
				ldx = ntop + nbot;
			}
		}
		if (ldx < ntop + nbot) {
			throw 0;
		}
				
		matrix_copy (ntop, ntop, a, x, lda, ldb);
		matrix_copy (nbot, ntop, a + ntop + n, x + ntop, lda, ldb);
		matrix_copy (ntop, nbot, a + (ntop + n) * lda, x + ntop * ldb, lda, ldb);
		matrix_copy (nbot, nbot, a + (ntop + n) * (lda + 1), x + ntop * (ldb + 1), lda, ldb);
		
		matrix_factorize (n, n, am, ipiv, info, lda);
		
		matrix_solve (n, am, ipiv, al, info, ntop, lda, lda);
		matrix_solve (n, am, ipiv, ar, info, nbot, lda, lda);
		
#ifdef _MPI
		matrix_matrix_multiply (ntop, ntop, n, -1.0, at, al, 1.0, x, lda, lda, ldb);
		matrix_matrix_multiply (nbot, ntop, n, -1.0, ab, al, 1.0, x + ntop, lda, lda, ldb);
		matrix_matrix_multiply (ntop, nbot, n, -1.0, at, ar, 1.0, x + ntop * ldx, lda, lda, ldb);
		matrix_matrix_multiply (nbot, nbot, n, -1.0, ab, ar, 1.0, x + ntop * (ldx + 1), lda, lda, ldb);

		if (id == 0) {
			std::vector <datatype> buffer (n2tot);
			std::vector <int> displs;
			if (id == 0) {
				displs.resize (np);
				for (int i = 1; i < np; ++i) {
					displs [i] = displs [i - 1] + ns2 [i - 1];
				}
			}
			MPI::COMM_WORLD.Gatherv (x, (ntop + nbot) * (ntop + nbot), MPI::DOUBLE, &buffer [0], &ns2 [0], &displs [0], MPI::DOUBLE, 0);

			for (int i = 0; i < ntot; ++i) {
				scale (ntot, 0.0, x + i * ldx);
			}
			int bcur = 0, xcur = 0;
			for (int i = 0; i < np - 1; ++i) {
				for (int j = 0; j < ns [i] + ns [i + 1]; ++j) {
					for (int k = 0; k < ns [i] + ns [i + 1]; ++k) {
						x [xcur + k * ldx + j] += buffer [bcur + k * (ns [i] + ns [i + 1]) + j];
					}
				}
				bcur += (ns [i] + ns [i + 1]) * (ns [i] + ns [i + 1]);
				xcur += ns [i] * (ldx + 1);
			}
			for (int j = 0; j < ns [np - 1]; ++j) {
				for (int k = 0; k < ns [np - 1]; ++k) {
					x [xcur + k * ldx + j] += buffer [bcur + k * ns [np - 1] + j];
				}
			}
			
			matrix_factorize (ntot, ntot, x, xipiv, info, ldx);
		} else {
			MPI::COMM_WORLD.Gatherv (x, (ntop + nbot) * (ntop + nbot), mpi_type (&typeid (datatype)), NULL, NULL, NULL, mpi_type (&typeid (datatype)), 0);
		}
#endif
		/*
			TODO Check that these work w/o MPI
		*/
	}
	
	template <class datatype>
	void p_block_matrix_solve (int id, int np, int n, int ntop, int nbot, datatype* a, int* ipiv, datatype* b, datatype *x, int *xipiv, int *ns, int *info, int nrhs = 1, int lda = -1, int ldx = -1, int ldb = -1) {
		std::vector <int> nsp (np);
		datatype* am = a + ntop * (lda + 1);
		datatype* at = a + ntop * lda;
		datatype* ab = a + ntop * lda + ntop + n;
		datatype* ar = a + (ntop + n) * lda + ntop;
		datatype* al = a + ntop;
		
		int ldy = 0;
		if (id == 0) {
			for (int i = 0; i < np - 1; ++i) {
				ldy += ns [i];
				nsp [i] = (ns [i] + ns [i + 1]) * nrhs;
			}
			ldy += ns [np - 1];
			nsp [np - 1] = (ns [np - 1]) * nrhs;
		} else {
			ldy = ntop + nbot;
		}
		
		std::vector <datatype> y (ldy * nrhs);
		
		if (lda == -1) {
			lda = n + ntop + nbot;
		}
		if (ldb == -1) {
			ldb = n + ntop + nbot;
		}
		if (ldx == -1) {
			ldx = ldy;
		}

		matrix_add_scaled (ntop, nrhs, 1.0, b, &y [0], ldb, ntop + nbot);
		matrix_add_scaled (nbot, nrhs, 1.0, b + ntop + n, &y [ntop], ldb, ntop + nbot);
		
		matrix_solve (n, am, ipiv, b + ntop, info, nrhs, lda, ldb);
				
#ifdef _MPI
		matrix_matrix_multiply (ntop, nrhs, n, -1.0, at, b + ntop, 1.0, &y [0], lda, ldb, ntop + nbot);
		matrix_matrix_multiply (nbot, nrhs, n, -1.0, ab, b + ntop, 1.0, &y [ntop], lda, ldb, ntop + nbot);
			
		if (id == 0) {
			int ycur = 0, bcur = 0;
			std::vector <datatype> buffer (nrhs * 2 * ldy);
			std::vector <int> displs;
			if (id == 0) {
				displs.resize (np);
				for (int i = 1; i < np; ++i) {
					displs [i] = displs [i - 1] + nsp [i - 1];
				}
			}
			MPI::COMM_WORLD.Gatherv (&y [0], (ntop + nbot) * nrhs, mpi_type (&typeid (datatype)), &buffer [0], &nsp [0], &displs [0], mpi_type (&typeid (datatype)), 0);
			utils::scale (nrhs * ldy, 0.0, &y [0]);
			for (int i = 0; i < np; ++i) {
				for (int j = 0; j < nrhs; ++j) {
					for (int k = 0; k < nsp [i] / nrhs; ++k) {
						y [j * ldy + ycur + k] += buffer [bcur + j * nsp [i] / nrhs + k];
					}
				}
				ycur += ns [i];
				bcur += nsp [i];
			}
			
			matrix_solve (ldy, x, xipiv, &y [0], info, nrhs, ldx, ldy);
			
			ycur = 0;
			bcur = 0;
			for (int i = 0; i < np; ++i) {
				for (int j = 0; j < nrhs; ++j) {
					for (int k = 0; k < nsp [i] / nrhs; ++k) {
						buffer [bcur + j * nsp [i] / nrhs + k] = y [j * ldy + ycur + k];
					}
				}
				ycur += ns [i];
				bcur += nsp [i];
			}
			
			MPI::COMM_WORLD.Scatterv (&buffer [0], &nsp [0], &displs [0], mpi_type (&typeid (datatype)), &y [0], (ntop + nbot) * nrhs, mpi_type (&typeid (datatype)), 0);
		} else {
			MPI::COMM_WORLD.Gatherv (&y [0], (ntop + nbot) * nrhs, mpi_type (&typeid (datatype)), NULL, NULL, NULL, mpi_type (&typeid (datatype)), 0);
			MPI::COMM_WORLD.Scatterv (NULL, NULL, NULL, mpi_type (&typeid (datatype)), &y [0], (ntop + nbot) * nrhs, mpi_type (&typeid (datatype)), 0);
		}
		
		matrix_copy (ntop, nrhs, &y [0], b, ntop + nbot, ldb);
		matrix_copy (nbot, nrhs, &y [ntop], b + ntop + n, ntop + nbot, ldb);
		
		matrix_matrix_multiply (n, nrhs, ntop, -1.0, al, b, 1.0, b + ntop, lda, ldb, ldb);
		matrix_matrix_multiply (n, nrhs, nbot, -1.0, ar, b + ntop + n, 1.0, b + ntop, lda, ldb, ldb);
#endif
	}
	
	void p_block_tridiag_factorize (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, int* ipiv, double *x, int *xipiv, int *info, int nrhs = 1, int lda = -1);
	
	void p_block_tridiag_solve (int id, int np, int n, double* sub, double *diag, double *sup, double *supsup, int* ipiv, double* b, double *x, int *xipiv, int *info, int nrhs = 1, int lda = -1, int ldb = -1);
	
	void p_block_banded_factorize (int id, int np, int n, int kl, int ku, double* matrix, int* ipiv, double *x, int *xipiv, double *bufferl, double *bufferr, int *info = NULL, int nrhs = 1, int lda = -1, int ldaa = -1);
	
	void p_block_banded_solve (int id, int np, int n, int kl, int ku, double* matrix, int* ipiv, double* b, double *x, int *xipiv, double *bufferl, double *bufferr, int *info = NULL, int nrhs = 1, int lda = -1, int ldaa = -1, int ldb = -1);
	
} /* utils */

#endif /* end of include guard: BLOCK_SOLVER_HPP_8C3ZNSDI */
