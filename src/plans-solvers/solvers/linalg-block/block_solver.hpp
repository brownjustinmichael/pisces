/*!**********************************************************************
 * \file block_solver.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_HPP_6489A25E
#define SOLVER_HPP_6489A25E

#include <vector>
#ifdef _MPI
#include <mpi.h>
#endif
#include "linalg/utils.hpp"
#include "linalg/linalg.hpp"
#include "mpi/messenger.hpp"

namespace linalg
{
	namespace block
	{
		/*!**********************************************************************
		 * \brief Factorize a block matrix
		 * 
		 * \param id The mpi id of this process
		 * \param np The number of mpi processes/number of blocks in the matrix
		 * \param n The number of elements in the non-overlapping part of the matrix
		 * \param ntop The top overlap region
		 * \param nbot The bottom overlap region
		 * \param a The matrix that holds the data (must be at least (n + ntop + nbot)^2)
		 * \param ipiv The integer array to hold the LAPACK swap information
		 * \param x Some additional memory needed for the solve, must be (ntop + nbot)^2 on all but process 0, for which it must be (sum (ntop + nbot))^2 over all processes
		 * \param xipiv The integer array to hold the LAPACK swap information for x
		 * \param ns An array of the integer ntop + nbot for all processes, only needed for process 0
		 * \param info A pointer to the result state of the solve, 0 if successful
		 * \param lda The leading dimension of a, if < 0, n + ntop + nbot
		 * \param ldx The leading dimension of x, if < 0, ntop + nbot for all but process 0
		 ************************************************************************/
		template <class datatype>
		void matrix_factorize (int id, int np, int n, int ntop, int nbot, datatype* a, int* ipiv, datatype *x, int *xipiv, int* ns, int *info, int lda = -1, int ldx = -1) {
			int ntot = 0, n2tot = 0, ldb = ntop + nbot;
			
			// Define some pointers for convenience
			datatype* am = a + ntop * (lda + 1);
			datatype* at = a + ntop * lda;
			datatype* ab = a + ntop * lda + ntop + n;
			datatype* ar = a + (ntop + n) * lda + ntop;
			datatype* al = a + ntop;
			
			if (lda == -1) {
				lda = n + ntop + nbot;
			}
			
			// Calculate the displacement array for the gather
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
				FATAL ("Leading dimension of x too small.");
				throw 0;
			}
			
			// Copy the corners of the matrix into x
			linalg::matrix_copy (ntop, ntop, a, x, lda, ldb);
			linalg::matrix_copy (nbot, ntop, a + ntop + n, x + ntop, lda, ldb);
			linalg::matrix_copy (ntop, nbot, a + (ntop + n) * lda, x + ntop * ldb, lda, ldb);
			linalg::matrix_copy (nbot, nbot, a + (ntop + n) * (lda + 1), x + ntop * (ldb + 1), lda, ldb);
			
			// Factorize the main matrix
			linalg::matrix_factorize (n, n, am, ipiv, info, lda);
			
			// Solve the matrix on the left and right sides
			linalg::matrix_solve (n, am, ipiv, al, info, ntop, lda, lda);
			linalg::matrix_solve (n, am, ipiv, ar, info, nbot, lda, lda);
			
	#ifdef _MPI
			// Multiply the left and right solutions by the top and bottom and add these into x
			linalg::matrix_matrix_multiply (ntop, ntop, n, -1.0, at, al, 1.0, x, lda, lda, ldb);
			linalg::matrix_matrix_multiply (nbot, ntop, n, -1.0, ab, al, 1.0, x + ntop, lda, lda, ldb);
			linalg::matrix_matrix_multiply (ntop, nbot, n, -1.0, at, ar, 1.0, x + ntop * ldx, lda, lda, ldb);
			linalg::matrix_matrix_multiply (nbot, nbot, n, -1.0, ab, ar, 1.0, x + ntop * (ldx + 1), lda, lda, ldb);

			if (id == 0) {
				// Construct the full displacement vector and gather x
				std::vector <datatype> buffer (n2tot);
				std::vector <int> displs;
				if (id == 0) {
					displs.resize (np);
					for (int i = 1; i < np; ++i) {
						displs [i] = displs [i - 1] + ns2 [i - 1];
					}
				}
				MPI::COMM_WORLD.Gatherv (x, (ntop + nbot) * (ntop + nbot), MPI::DOUBLE, &buffer [0], &ns2 [0], &displs [0], MPI::DOUBLE, 0);
				
				// Reconstitute the matrix (some of the elements overlap and need to be summed)
				for (int i = 0; i < ntot; ++i) {
					linalg::scale (ntot, 0.0, x + i * ldx);
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
				
				// Factorize the boundary matrix
				linalg::matrix_factorize (ntot, ntot, x, xipiv, info, ldx);
			} else {
				// Send x to process 0 and exit
				MPI::COMM_WORLD.Gatherv (x, (ntop + nbot) * (ntop + nbot), mpi::mpi_type (&typeid (datatype)), NULL, NULL, NULL, mpi::mpi_type (&typeid (datatype)), 0);
			}
	#endif
			/*
				TODO Check that these work w/o MPI
			*/
		}
		
		/*!**********************************************************************
		 * \brief Solve a block matrix
		 * 
		 * \param id The mpi id of this process
		 * \param np The number of mpi processes/number of blocks in the matrix
		 * \param n The number of elements in the non-overlapping part of the matrix
		 * \param ntop The top overlap region
		 * \param nbot The bottom overlap region
		 * \param a The matrix that holds the data (must be at least (n + ntop + nbot)^2)
		 * \param ipiv The integer array to hold the LAPACK swap information
		 * \param b The right hand side to solve (overwritten with solution)
		 * \param x Some additional memory needed for the solve, must be (ntop + nbot)^2 on all but process 0, for which it must be (sum (ntop + nbot))^2 over all processes
		 * \param xipiv The integer array to hold the LAPACK swap information for x
		 * \param ns An array of the integer ntop + nbot for all processes, only needed for process 0
		 * \param info A pointer to the result state of the solve, 0 if successful
		 * \param nrhs The number of right hand sides in b
		 * \param lda The leading dimension of a, if < 0, n + ntop + nbot
		 * \param ldx The leading dimension of x, if < 0, ntop + nbot for all but process 0
		 * \param ldb The leading dimension of b, if < 0, n + ntop + nbot
		 ************************************************************************/
		template <class datatype>
		void matrix_solve (int id, int np, int n, int ntop, int nbot, datatype* a, int* ipiv, datatype* b, datatype *x, int *xipiv, int *ns, int *info, int nrhs = 1, int lda = -1, int ldx = -1, int ldb = -1) {
			std::vector <int> nsp (np);
			
			// Define some pointers in a for convenience
			datatype* am = a + ntop * (lda + 1);
			datatype* at = a + ntop * lda;
			datatype* ab = a + ntop * lda + ntop + n;
			datatype* ar = a + (ntop + n) * lda + ntop;
			datatype* al = a + ntop;
			
			// We'll need a right hand side for x, which we'll call y. Calculate its dimensions
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
			
			// Copy the top and bottom of b into y
			linalg::matrix_add_scaled (ntop, nrhs, 1.0, b, &y [0], ldb, ntop + nbot);
			linalg::matrix_add_scaled (nbot, nrhs, 1.0, b + ntop + n, &y [ntop], ldb, ntop + nbot);
			
			// Solve the main body of b
			linalg::matrix_solve (n, am, ipiv, b + ntop, info, nrhs, lda, ldb);
				
	#ifdef _MPI
			// Multiply b by the top and bottom of a
			linalg::matrix_matrix_multiply (ntop, nrhs, n, -1.0, at, b + ntop, 1.0, &y [0], lda, ldb, ntop + nbot);
			linalg::matrix_matrix_multiply (nbot, nrhs, n, -1.0, ab, b + ntop, 1.0, &y [ntop], lda, ldb, ntop + nbot);
			
			if (id == 0) {
				int ycur = 0, bcur = 0;
				std::vector <datatype> buffer (nrhs * 2 * ldy);
				std::vector <int> displs;
				// Calculate the true displacement vector for the gather
				if (id == 0) {
					displs.resize (np);
					for (int i = 1; i < np; ++i) {
						displs [i] = displs [i - 1] + nsp [i - 1];
					}
				}
				MPI::COMM_WORLD.Gatherv (&y [0], (ntop + nbot) * nrhs, mpi::mpi_type (&typeid (datatype)), &buffer [0], &nsp [0], &displs [0], mpi::mpi_type (&typeid (datatype)), 0);
				linalg::scale (nrhs * ldy, 0.0, &y [0]);
				// Some of the values will overlap, as before, so we compose y
				for (int i = 0; i < np; ++i) {
					for (int j = 0; j < nrhs; ++j) {
						for (int k = 0; k < nsp [i] / nrhs; ++k) {
							y [j * ldy + ycur + k] += buffer [bcur + j * nsp [i] / nrhs + k];
						}
					}
					ycur += ns [i];
					bcur += nsp [i];
				}
				
				// Solve y
				linalg::matrix_solve (ldy, x, xipiv, &y [0], info, nrhs, ldx, ldy);
			
				ycur = 0;
				bcur = 0;
				// We reconstitute y to send back to the other processes
				for (int i = 0; i < np; ++i) {
					for (int j = 0; j < nrhs; ++j) {
						for (int k = 0; k < nsp [i] / nrhs; ++k) {
							buffer [bcur + j * nsp [i] / nrhs + k] = y [j * ldy + ycur + k];
						}
					}
					ycur += ns [i];
					bcur += nsp [i];
				}
				
				MPI::COMM_WORLD.Scatterv (&buffer [0], &nsp [0], &displs [0], mpi::mpi_type (&typeid (datatype)), &y [0], (ntop + nbot) * nrhs, mpi::mpi_type (&typeid (datatype)), 0);
			} else {
				// Send the edge values and wait to receive the solved edge values
				MPI::COMM_WORLD.Gatherv (&y [0], (ntop + nbot) * nrhs, mpi::mpi_type (&typeid (datatype)), NULL, NULL, NULL, mpi::mpi_type (&typeid (datatype)), 0);
				MPI::COMM_WORLD.Scatterv (NULL, NULL, NULL, mpi::mpi_type (&typeid (datatype)), &y [0], (ntop + nbot) * nrhs, mpi::mpi_type (&typeid (datatype)), 0);
			}
			
			// Copy the received values into b
			linalg::matrix_copy (ntop, nrhs, &y [0], b, ntop + nbot, ldb);
			linalg::matrix_copy (nbot, nrhs, &y [ntop], b + ntop + n, ntop + nbot, ldb);
			
			// The values in the body of b will need to be updated now that the edge values are known
			linalg::matrix_matrix_multiply (n, nrhs, ntop, -1.0, al, b, 1.0, b + ntop, lda, ldb, ldb);
			linalg::matrix_matrix_multiply (n, nrhs, nbot, -1.0, ar, b + ntop + n, 1.0, b + ntop, lda, ldb, ldb);
	#endif
		}
	} /* block */
} /* linalg */

#endif /* end of include guard: SOLVER_HPP_6489A25E */
