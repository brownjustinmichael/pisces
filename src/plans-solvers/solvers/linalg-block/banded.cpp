/*!**********************************************************************
 * \file block_solver.cpp
 * /Users/justinbrown/Dropbox/pisces/src/utils
 * 
 * Created by Justin Brown on 2013-11-17.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <vector>
#include <sstream>
#ifdef _MPI
#include <mpi.h>
#endif
#include "mpi/messenger.hpp"
#include "linalg/utils.hpp"
#include "linalg/linalg.hpp"
#include "banded.hpp"

namespace linalg
{
	void p_block_banded_factorize (int id, int np, int n, int kl, int ku, double* matrix, int* ipiv, double *x, int *xipiv, double *bufferl, double *bufferr, int *info, int nrhs, int lda, int ldaa) {
		int ntop, nbot;
		std::stringstream debug;
		int ldx = 2 * (kl + ku);
		int ldxx = ldx * ldx;
		if (id == 0) {
			ntop = 0;
			ldxx *= np * np;
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
		if (ldaa == -1) {
			ldaa = n + ku + kl + ntop + nbot;
		}
		
		linalg::scale (nrhs * ldxx, 0.0, x);
		linalg::scale (nrhs * kl * n, 0.0, bufferl);
		linalg::scale (nrhs * ku * n, 0.0, bufferr);
		
		for (int i = 0; i < nrhs; ++i) {
			// for (int j = kl + ntop; j < n + kl + ntop; ++j) {
			// 	for (int k = 0; k < lda; ++k) {
			// 		debug << matrix [i * lda * ldaa + j * lda + k] << " ";
			// 	}
			// 	DEBUG ("MAT " << debug.str ());
			// 	debug.str ("");
			// }
			linalg::matrix_banded_factorize (n, n, kl, ku, matrix + (i) * lda * ldaa + (kl + ntop) * lda, ipiv + i * n, info, lda);
		}

#ifdef _MPI
		if (id != 0) {
			for (int i = 0; i < nrhs; ++i) {
				for (int j = 0; j < kl + ku; ++j) {
					for (int k = 0; k < std::min (ku, j + 1); ++k) {
						x [i * ldxx + j * ldx + k + kl] = matrix [i * ldaa * lda + j * lda + 2 * kl + ku + k - j];
					}
				}
				for (int j = 0; j < kl; ++j) {
					for (int k = 0; k < j + 1; ++k) {
						bufferl [i * n * kl + j * n + k] = matrix [i * lda * ldaa + (j + ku) * lda + k + 2 * kl + ku - j];
					}
				}

				linalg::matrix_banded_solve (n, kl, ku, matrix + i * lda * ldaa + (ntop + kl) * lda, ipiv + i * n, &bufferl [i * n * kl], info, kl, lda);

				for (int j = 0; j < kl; ++j) {
					for (int k = 0; k < ku; ++k) {
						for (int l = 0; l < k + 1; ++l) {
							x [i * ldxx + (j + ku) * ldx + k + kl] -= bufferl [i * n * kl + j * n + l + k] * matrix [i * lda * ldaa + (k + kl + ntop + l) * lda + kl - l];
						}
					}
					if (id != np - 1) {
						for (int k = 0; k < kl; ++k) {
							for (int l = 0; l < kl - k; ++l) {
								x [i * ldxx + (j + ku) * ldx + k + kl + ku] -= bufferl [i * n * kl + j * n + l + k + n - kl] * matrix [i * ldaa * lda + (k + l + n + ntop) * lda + 2 * kl + ku - l];
							}
						}
					}
				}
			}
		}

		if (id != np - 1) {
			for (int i = 0; i < nrhs; ++i) {
				for (int j = 0; j < ku + kl; ++j) {
					for (int k = std::max (0, j - ku); k < kl; ++k) {
						x [i * ldxx + (j + ku + kl) * ldx + ku + kl + k] = matrix [i * ldaa * lda + (j + n + ntop + kl) * lda + kl + ku + k - j];
					}
				}
				for (int j = 0; j < ku; ++j) {
					for (int k = -1; k >= -ku; --k) {
						bufferr [i * n * ku + j * n + k + n] = matrix [i * lda * ldaa + (n + ntop + j + kl) * lda + k + kl + ku - j];
					}
				}

				linalg::matrix_banded_solve (n, kl, ku, matrix + i * lda * ldaa + (ntop + kl) * lda, ipiv + i * n, &bufferr [i * n * ku], info, ku, lda);

				for (int j = 0; j < ku; ++j) {
					if (id != 0) {
						for (int k = 0; k < ku; ++k) {
							for (int l = 0; l < k + 1; ++l) {
								x [i * ldxx + (j + ku + kl) * ldx + k + kl] -= bufferr [i * n * ku + j * n + l + k] * matrix [i * lda * ldaa + (k + kl + ntop + l) * lda + kl - l];
							}
						}
					}
					for (int k = 0; k < kl; ++k) {
						for (int l = 0; l < kl - k; ++l) {
							x [i * ldxx + (j + ku + kl) * ldx + k + kl + ku] -= bufferr [i * n * ku + j * n + l + k + n - kl] * matrix [i * ldaa * lda + (k + l + n + ntop) * lda + 2 * kl + ku - l];
						}
					}
				}
			}
		}

		if (id == 0) {
			ldx = np * (2 * ku + 2 * kl);
			std::vector <double> buffer (np * 4 * (ku + kl) * (ku + kl));
			for (int q = 0; q < nrhs; ++q) {
				MPI::COMM_WORLD.Gather (x + q * ldxx, 4 * (ku + kl) * (ku + kl), MPI::DOUBLE, &buffer [0], 4 * (ku + kl) * (ku + kl), MPI::DOUBLE, 0);

				linalg::scale (ldx * ldx, 0.0, x + q * ldxx);

				int bcur = 0, xcur = 0;
				for (int i = 0; i < np - 1; ++i) {
					for (int j = 0; j < 2 * ku + 2 * kl; ++j) {
						for (int k = 0; k < 2 * kl + 2 * ku; ++k) {
							x [q * ldxx + xcur + k * ldx + j] += buffer [bcur + k * 2 * (ku + kl) + j];
						}
					}
					bcur += (2 * ku + 2 * kl) * 2 * (ku + kl);
					xcur += (ku + kl) * (ldx + 1);
				}
				for (int j = 0; j < 2 * ku + 2 * kl; ++j) {
					for (int k = 0; k < 2 * kl + 2 * ku; ++k) {
						x [q * ldxx + xcur + k * ldx + j] += buffer [bcur + k * 2 * (ku + kl) + j];
					}
				}
				
				linalg::matrix_factorize ((ku + kl) * (np - 1), (ku + kl) * (np - 1), x + q * ldxx + (kl + ku) * (ldx + 1), xipiv + q * ldx, info, ldx);
			}

		} else {
			for (int q = 0; q < nrhs; ++q) {
				MPI::COMM_WORLD.Gather (x + q * ldxx, (2 * ku + 2 * kl) * 2 * (ku + kl), mpi::mpi_type (&typeid (double)), NULL, (2 * ku + 2 * kl) * 2 * (ku + kl), mpi::mpi_type (&typeid (double)), 0);
			}
		}
#endif
	}

	void p_block_banded_solve (int id, int np, int n, int kl, int ku, double* matrix, int* ipiv, double* b, double *x, int *xipiv, double *bufferl, double *bufferr, int *info, int nrhs, int lda, int ldaa, int ldb) {
		std::stringstream debug;
		std::vector <double> y (2 * (ku + kl) * np * nrhs, 0.0);
		int ntop, nbot;
		int ldx = 2 * (kl + ku);
		int ldxx = ldx * ldx;
		if (id == 0) {
			ntop = 0;
			ldxx *= np * np;
		} else {
			ntop = ku;
		}
		if (id == np - 1) {
			nbot = 0;
		} else {
			nbot = kl;
		}
		
		int ldy;
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
		if (ldaa == -1) {
			ldaa = n + ku + kl + ntop + nbot;
		}

		linalg::matrix_add_scaled (ntop, nrhs, 1.0, b, &y [kl], ldb, ldy);
		linalg::matrix_add_scaled (nbot, nrhs, 1.0, b + ntop + n, &y [kl + ku], ldb, ldy);
		
		for (int i = 0; i < nrhs; ++i) {
			linalg::matrix_banded_solve (n, kl, ku, matrix + i * lda * ldaa + (kl + ntop) * lda, ipiv + i * n, b + ntop + i * ldb, info, 1, lda);
		}
		
#ifdef _MPI
		if (id != 0) {
			for (int i = 0; i < nrhs; ++i) {
				for (int k = 0; k < ku; ++k) {
					for (int l = 0; l < k + 1; ++l) {
						y [i * ldy + k + kl] -= b [i * ldb + l + k + ntop] * matrix [i * lda * ldaa + (k + kl + ntop + l) * lda + kl - l];
					}
				}
			}
		}

		if (id != np - 1) {
			for (int i = 0; i < nrhs; ++i) {
				for (int k = 0; k < kl; ++k) {
					for (int l = 0; l < kl - k; ++l) {
						y [i * ldy + k + kl + ku] -= b [i * ldb + l + k + ntop + n - kl] * matrix [i * ldaa * lda + (k + l + n + ntop) * lda + 2 * kl + ku - l];
					}
				}
			}
		}
		if (id == 0) {
			int ycur = 0, bcur = 0;
			ldx = np * (2 * ku + 2 * kl);

			std::vector <double> buffer (2 * (ku + kl) * np);
			for (int j = 0; j < nrhs; ++j) {
				ycur = 0;
				bcur = 0;
				
				MPI::COMM_WORLD.Gather (&y [j * ldy], 2 * (kl + ku), mpi::mpi_type (&typeid (double)), &buffer [0], 2 * (kl + ku), mpi::mpi_type (&typeid (double)), 0);
				linalg::scale (ldy, 0.0, &y [j * ldy]);
				for (int i = 0; i < np; ++i) {
						for (int k = 0; k < 2 * (ku + kl); ++k) {
							y [j * ldy + ycur + k] += buffer [bcur + k];
						}
					ycur += ku + kl;
					bcur += 2 * (ku + kl);
				}

				linalg::matrix_solve ((ku + kl) * (np - 1), x + j * ldx * ldx + (kl + ku) * (ldx + 1), xipiv + j * ldx, &y [j * ldy + kl + ku], info, 1, ldx);

				ycur = 0;
				bcur = 0;
				for (int i = 0; i < np; ++i) {
					for (int k = 0; k < 2 * (ku + kl); ++k) {
						buffer [bcur + k] = y [j * ldy + ycur + k];
					}
					ycur += ku + kl;
					bcur += 2 * (ku + kl);
				}
				
				MPI::COMM_WORLD.Scatter (&buffer [0], 2 * (kl + ku), mpi::mpi_type (&typeid (double)), &y [j * ldy], 2 * (kl + ku), mpi::mpi_type (&typeid (double)), 0);
			}
		} else {
			for (int j = 0; j < nrhs; ++j) {
				MPI::COMM_WORLD.Gather (&y [j * ldy], 2 * (kl + ku), mpi::mpi_type (&typeid (double)), NULL, 2 * (kl + ku), mpi::mpi_type (&typeid (double)), 0);
				MPI::COMM_WORLD.Scatter (NULL, 2 * (kl + ku), mpi::mpi_type (&typeid (double)), &y [j * ldy], 2 * (kl + ku), mpi::mpi_type (&typeid (double)), 0);
			}
		}
		
		if (id != 0) {
			linalg::matrix_copy (ntop, nrhs, &y [kl], b, ldy, ldb);
		}
		if (id != np - 1) {
			linalg::matrix_copy (nbot, nrhs, &y [kl + ku], b + ntop + n, ldy, ldb);
		}
		
		for (int i = 0; i < nrhs; ++i) {
			if (id != 0) {
				linalg::matrix_matrix_multiply (n, 1, kl, -1.0, &bufferl [i * n * kl], &y [ku + i * ldy], 1.0, b + ntop + i * ldb, n, ldb, ldb);
			}
			if (id != np - 1) {
				linalg::matrix_matrix_multiply (n, 1, ku, -1.0, &bufferr [i * n * ku], &y [kl + ku + i * ldy], 1.0, b + ntop + i * ldb, n, ldb, ldb);
			}
		}
#endif
	}
} /* linalg */