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
#include "linalg/utils.hpp"
#include "linalg/linalg.hpp"
#include "block_solver.hpp"

namespace utils
{
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
		
		utils::scale (nrhs * ldxx, 0.0, x);
		utils::scale (nrhs * kl * n, 0.0, bufferl);
		utils::scale (nrhs * ku * n, 0.0, bufferr);
		
		DEBUG (n << " " << lda << " " << ldaa)

		for (int i = 0; i < nrhs; ++i) {
			for (int k = 0; k < 2 * kl + ku + 1; ++k) {
				for (int j = 0; j < n; ++j) {
					debug << matrix [i * lda * ldaa + (kl + ntop + j) * lda + k] << " ";
				}
				DEBUG (debug.str ());
				debug.str ("");
			}
			matrix_banded_factorize (n, n, kl, ku, matrix + (i) * lda * ldaa + (kl + ntop) * lda, ipiv + i * n, info, lda);
			for (int k = 0; k < 2 * kl + ku + 1; ++k) {
				for (int j = 0; j < n; ++j) {
					debug << matrix [i * lda * ldaa + (kl + ntop + j) * lda + k] << " ";
				}
				DEBUG ("AFTER " << debug.str ());
				debug.str ("");
			}
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

				matrix_banded_solve (n, kl, ku, matrix + i * lda * ldaa + (ntop + kl) * lda, ipiv + i * n, &bufferl [i * n * kl], info, kl, lda);

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

				matrix_banded_solve (n, kl, ku, matrix + i * lda * ldaa + (ntop + kl) * lda, ipiv + i * n, &bufferr [i * n * ku], info, ku, lda);

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

				scale (ldx * ldx, 0.0, x + q * ldxx);

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
				
				for (int j = 0; j < (ku + kl) * (np - 1); ++j) {
					for (int k = 0; k < (ku + kl) * (np - 1); ++k) {
						debug << x [q * ldxx + (kl + ku) * (ldx + 1) + j * ldx + k] << " ";
					}
					DEBUG (debug.str ());
					debug.str ("");
				}
				
				matrix_factorize ((ku + kl) * (np - 1), (ku + kl) * (np - 1), x + q * ldxx + (kl + ku) * (ldx + 1), xipiv + q * ldx, info, ldx);
			}

		} else {
			for (int q = 0; q < nrhs; ++q) {
				MPI::COMM_WORLD.Gather (x + q * ldxx, (2 * ku + 2 * kl) * 2 * (ku + kl), mpi_type (&typeid (double)), NULL, (2 * ku + 2 * kl) * 2 * (ku + kl), mpi_type (&typeid (double)), 0);
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

		matrix_add_scaled (ntop, nrhs, 1.0, b, &y [kl], ldb, ldy);
		matrix_add_scaled (nbot, nrhs, 1.0, b + ntop + n, &y [kl + ku], ldb, ldy);
		
		DEBUG (n << " " << lda << " " << ldaa);

		for (int i = 0; i < nrhs; ++i) {
			for (int j = 0; j < n + ntop + nbot; ++j) {
				debug << b [i * ldb + j] << " ";
			}
			DEBUG (debug.str ());
			debug.str ("");
			for (int k = 0; k < 2 * kl + ku + 1; ++k) {
				for (int j = 0; j < n; ++j) {
					debug << matrix [i * lda * ldaa + (kl + ntop + j) * lda + k] << " ";
				}
				DEBUG ("MATRIX " << debug.str ());
				debug.str ("");
			}
			matrix_banded_solve (n, kl, ku, matrix + i * lda * ldaa + (kl + ntop) * lda, ipiv + i * n, b + ntop + i * ldb, info, 1, lda);
			for (int j = 0; j < n + ntop + nbot; ++j) {
				debug << b [i * ldb + j] << " ";
			}
			DEBUG ("BEFORE " << debug.str ());
			debug.str ("");
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
				
				MPI::COMM_WORLD.Gather (&y [j * ldy], 2 * (kl + ku), mpi_type (&typeid (double)), &buffer [0], 2 * (kl + ku), mpi_type (&typeid (double)), 0);
				utils::scale (ldy, 0.0, &y [j * ldy]);
				for (int i = 0; i < np; ++i) {
						for (int k = 0; k < 2 * (ku + kl); ++k) {
							y [j * ldy + ycur + k] += buffer [bcur + k];
						}
					ycur += ku + kl;
					bcur += 2 * (ku + kl);
				}

				matrix_solve ((ku + kl) * (np - 1), x + j * ldx * ldx + (kl + ku) * (ldx + 1), xipiv + j * ldx, &y [j * ldy + kl + ku], info, 1, ldx);

				ycur = 0;
				bcur = 0;
				for (int i = 0; i < np; ++i) {
					for (int k = 0; k < 2 * (ku + kl); ++k) {
						buffer [bcur + k] = y [j * ldy + ycur + k];
					}
					ycur += ku + kl;
					bcur += 2 * (ku + kl);
				}
				
				MPI::COMM_WORLD.Scatter (&buffer [0], 2 * (kl + ku), mpi_type (&typeid (double)), &y [j * ldy], 2 * (kl + ku), mpi_type (&typeid (double)), 0);
			}
		} else {
			for (int j = 0; j < nrhs; ++j) {
				MPI::COMM_WORLD.Gather (&y [j * ldy], 2 * (kl + ku), mpi_type (&typeid (double)), NULL, 2 * (kl + ku), mpi_type (&typeid (double)), 0);
				MPI::COMM_WORLD.Scatter (NULL, 2 * (kl + ku), mpi_type (&typeid (double)), &y [j * ldy], 2 * (kl + ku), mpi_type (&typeid (double)), 0);
			}
		}
		
		if (id != 0) {
			matrix_copy (ntop, nrhs, &y [kl], b, ldy, ldb);
		}
		if (id != np - 1) {
			matrix_copy (nbot, nrhs, &y [kl + ku], b + ntop + n, ldy, ldb);
		}
		
		for (int i = 0; i < nrhs; ++i) {
			if (id != 0) {
				matrix_matrix_multiply (n, 1, kl, -1.0, &bufferl [i * n * kl], &y [ku + i * ldy], 1.0, b + ntop + i * ldb, n, ldb, ldb);
			}
			if (id != np - 1) {
				matrix_matrix_multiply (n, 1, ku, -1.0, &bufferr [i * n * ku], &y [kl + ku + i * ldy], 1.0, b + ntop + i * ldb, n, ldb, ldb);
			}
		}
#endif
	}
} /* utils */