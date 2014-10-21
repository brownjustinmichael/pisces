/*!**********************************************************************
 * \file block_test.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-11-05.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cxxtest/TestSuite.h>

#include "mpi/messenger.hpp"
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include "linalg/linalg.hpp"
#include "plans-solvers/solvers/linalg-block/tridiagonal.hpp"
#include "plans-solvers/solvers/linalg-block/banded.hpp"
#include "plans-solvers/solvers/linalg-block/solver.hpp"
#include <iostream>

#define TEST_TINY 1.e-6

class mpi_solver_test_suite : public CxxTest::TestSuite
{
private:
	mpi::messenger mess;
	
public:
	void test_block_solver () {
		int timesteps = 1;
		int id = mess.get_id ();
		int np = mess.get_np ();
		int n = 10, nrhs = 10, ntop = 0, nbot = 0;
		if (id != 0) {
			ntop = 1;
		}
		if (id != np - 1) {
			nbot = 1;
		}
		int lda = n + ntop + nbot, ldx = 0, ldb = lda;
	
		std::vector <int> ns (mess.get_np ());
		mess.gather <int> (1, &ntop, &ns [0]);
		if (mess.get_id () == 0) {
			for (int i = 0; i < mess.get_np (); ++i) {
				ldx += ns [i];
			}
		} else {
			ldx = ntop + nbot;
		}
		
		std::vector <double> bufftop (ntop * nrhs), buffbot (nbot * nrhs), rbufftop (ntop * nrhs), rbuffbot (nbot * nrhs);
	
		std::vector <double> a (lda * lda), acopy (lda * lda);
		std::vector <double> b (lda * nrhs), bcopy (lda * nrhs);
		std::vector <int> ipiv (lda);
		std::vector <double> x (ldx * ldx);
		std::vector <int> xipiv (ldx);
		int info;
	
		srand (2);
	
		for (int i = 0; i < lda; ++i) {
			for (int j = 0; j < lda; ++j) {
				a [j * lda + i] = rand () % 100;
			}
			for (int j = 0; j < nrhs; ++j) {
				b [j * ldb + i] = rand () % 100;
				bcopy [j * ldb + i] = b [j * ldb + i];
			}
		}
	
		linalg::matrix_copy (lda, lda, &a [0], &acopy [0], lda, lda);
	
		clock_t corig_begin, corig_end;
		std::chrono::time_point <std::chrono::system_clock> orig_begin, orig_end;
	
		corig_begin = clock ();
		orig_begin = std::chrono::system_clock::now ();
	
		linalg::matrix_factorize (n, n, &a [0], &ipiv [0], &info, lda);
		linalg::matrix_solve (n, &a [0], &ipiv [0], &b [0], &info, nrhs, lda, ldb);
	
		corig_end = clock ();
		orig_end = std::chrono::system_clock::now ();
	
		std::chrono::duration <double> orig_total = orig_end - orig_begin;
	
		// printf ("[%i] Orig CPU Time: %f\n", id, ((double) (corig_end - corig_begin))/CLOCKS_PER_SEC);
		// printf ("[%i] Orig Wall Time: %f\n", id, (double) orig_total.count ());
	
		clock_t cbegin, cmid, cend;
		std::chrono::time_point <std::chrono::system_clock> begin, mid, end;
	
		std::chrono::duration <double> mb, em, eb;
		double cmb = 0.0, cem = 0.0, ceb = 0.0;
	
		for (int i = 0; i < timesteps; ++i) {
			for (int i = 0; i < lda; ++i) {
				for (int j = 0; j < lda; ++j) {
					a [j * lda + i] = rand () % 100;
					acopy [j * lda + i] = a [j * lda + i];
				}
				for (int j = 0; j < nrhs; ++j) {
					b [j * ldb + i] = rand () % 100;
					bcopy [j * ldb + i] = b [j * ldb + i];
				}
			}
			
			cbegin = clock ();
			begin = std::chrono::system_clock::now ();
		
			linalg::p_block_matrix_factorize (mess.get_id (), mess.get_np (), n, ntop, nbot, &a [0], &ipiv [0], &x [0], &xipiv [0], &ns [0], &info, lda, ldx);
	
			cmid = clock ();
			mid = std::chrono::system_clock::now ();
		
			mb += mid - begin;
			cmb += cmid - cbegin;
		
			linalg::p_block_matrix_solve (mess.get_id (), mess.get_np (), n, ntop, nbot, &a [0], &ipiv [0], &b [0], &x [0], &xipiv [0], &ns [0], &info, nrhs, lda, ldx, ldb);
	
			cend = clock ();
			end = std::chrono::system_clock::now ();
		
			em += end - mid;
			eb += end - begin;
			cem += cend - cmid;
			ceb += cend - cbegin;
			
			linalg::matrix_matrix_multiply (n + ntop + nbot, nrhs, n + ntop + nbot, 1.0, &acopy [0], &b [0], -1.0, &bcopy [0]);
			
			if (id != 0) {
				linalg::matrix_copy (ntop, nrhs, &bcopy [0], &bufftop [0], ldb, ntop);
				mess.send (ntop * nrhs, &bufftop [0], id - 1, 0);
			}
			if (id != np - 1) {
				linalg::matrix_copy (nbot, nrhs, &bcopy [n + ntop], &buffbot [0], ldb, nbot);
				mess.recv (nbot * nrhs, &rbuffbot [0], id + 1, 0);
				mess.send (nbot * nrhs, &buffbot [0], id + 1, 1);
				linalg::matrix_add_scaled (nbot, nrhs, 1.0, &rbuffbot [0], &bcopy [n + ntop], nbot, ldb);
			}
			if (id != 0) {
				mess.recv (ntop * nrhs, &rbufftop [0], id - 1, 1);
				linalg::matrix_add_scaled (ntop, nrhs, 1.0, &rbufftop [0], &bcopy [0], ntop, ldb);
			}
			
			for (int i = 0; i < n + ntop + nbot; ++i) {
				for (int j = 0; j < nrhs; ++j) {
					TSM_ASSERT_DELTA ("Block solver failure", bcopy [j * ldb + i], 0.0, TEST_TINY);
				}
			}
		}
	
		// printf ("[%i] Avg CPU Time: %f + %f = %f\n", id, ((double) (cmb))/CLOCKS_PER_SEC/timesteps, ((double) (cem))/CLOCKS_PER_SEC/timesteps, ((double) (ceb))/CLOCKS_PER_SEC/timesteps);
		// printf ("[%i] Avg Wall Time: %f + %f = %f\n", id, (double) mb.count ()/timesteps, (double) em.count ()/timesteps, (double) eb.count ()/timesteps);
	}
	
	void test_block_tridiagonal_solve () {
		// io::parameters <double> params ("../input/block_parameters.txt");
		// omp_set_num_threads(params.nmp);

		int id = mess.get_id ();
		int np = mess.get_np ();
		int n = 60, nrhs = 70, ntop = 0, nbot = 0;
		if (id != 0) {
			ntop = 1;
		}
		if (id != np - 1) {
			nbot = 1;
		}

		std::vector <double> sub (n * nrhs, 0.0);
		std::vector <double> diag (n * nrhs, 0.0);
		std::vector <double> sup (n * nrhs, 0.0);
		std::vector <double> subcopy (n * nrhs, 0.0);
		std::vector <double> diagcopy (n * nrhs, 0.0);
		std::vector <double> supcopy (n * nrhs, 0.0);
		std::vector <double> supsup (n * nrhs, 0.0);
		std::vector <int> ipiv (n * nrhs, 0);

		std::vector <double> x (2 * n * nrhs + 8 * np * nrhs, 0.0);
		std::vector <int> xipiv (2 * np * nrhs, 0);

		std::vector <double> b (n * nrhs, 0.0);
		std::vector <double> bcopy (n * nrhs, 0.0);

		int info;

		srand (1);
		for (int i = 0; i < nrhs; ++i) {
			for (int j = 0; j < n; ++j) {
				sub [i * n + j] = rand () % 100;
				subcopy [i * n + j] = sub [i * n + j];
				diag [i * n + j] = rand () % 100;
				diagcopy [i * n + j] = diag [i * n + j];
				sup [i * n + j] = rand () % 100;
				supcopy [i * n + j] = sup [i * n + j];
				b [i * n + j] = rand () % 100;
				bcopy [i * n + j] = b [i * n + j];
			}
		}
		
		try {
			linalg::p_block_tridiag_factorize (id, np, n - ntop - nbot, &sub [0], &diag [0], &sup [0], &supsup [0], &ipiv [0], &x [0], &xipiv [0], &info, nrhs);

			linalg::p_block_tridiag_solve (id, np, n - ntop - nbot,  &sub [0], &diag [0], &sup [0], &supsup [0], &ipiv [0], &b [0], &x [0], &xipiv [0], &info, nrhs);
		} catch (std::exception& except) {
			std::cout << except.what () << '\n';
		}

		for (int j = 0; j < nrhs; ++j) {
			bcopy [j * n] -= diagcopy [j * n] * b [j * n] + supcopy [j * n] * b [1 + j * n];
			for (int i = 1; i < n - 1; ++i) {
				bcopy [i + j * n] -= subcopy [i + j * n] * b [i - 1 + j * n] + diagcopy [i + j * n] * b [i + j * n] + supcopy [i + j * n] * b [i + 1 + j * n];
			}
			bcopy [n - 1 + j * n] -= subcopy [n - 1 + j * n] * b [n - 2 + j * n] + diagcopy [n - 1 + j * n] * b [n - 1 + j * n];
		}
	
		std::vector <double> edge_0 (nrhs), edge_n (nrhs), redge_0 (nrhs), redge_n (nrhs);
	
		for (int i = 0; i < nrhs; ++i) {
			edge_0 [i] = b [i * n];
			edge_n [i] = b [i * n + n - 1];
		}

		if (id != 0) {
			mess.send (nrhs, &edge_0 [0], id - 1, 0);
		}
		if (id != np - 1) {
			mess.recv (nrhs, &redge_n [0], id + 1, 0);
			for (int i = 0; i < nrhs; ++i) {
				bcopy [i * n + n - 1] -= supcopy [i * n + n - 1] * redge_n [i];
			}
		}
		if (id != np - 1) {
			mess.send (nrhs, &edge_n [0], id + 1, 1);
		}
		if (id != 0) {
			mess.recv (nrhs, &redge_0 [0], id - 1, 1);
			for (int i = 0; i < nrhs; ++i) {
				bcopy [i * n] -= subcopy [i * n] * redge_0 [i];
			}
		}
		
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < nrhs; ++j) {
				TSM_ASSERT_DELTA ("Tridiagonal block solver failure", bcopy [j * n + i], 0.0, TEST_TINY);
			}
		}
	}
	
	void test_block_banded_solve () {
		int id = mess.get_id ();
		int np = mess.get_np ();
		int n = 60, nrhs = 10, ntop = 0, nbot = 0;

		int ku = 1, kl = 2;
		if (id != 0) {
			ntop = ku;
		}
		if (id != np - 1) {
			nbot = kl;
		}
		int lda = ku + 2 * kl + 1 + 4;
		int ldaa = n + ku + kl + 8;
		int ldb = ldaa + 10;

		std::vector <double> matrix (lda * ldaa * nrhs, 0.0);
		std::vector <double> matrixcopy (lda * ldaa * nrhs, 0.0);
		std::vector <double> bufferl (kl * nrhs * n, 0.0), bufferr (ku * nrhs * n, 0.0);
		std::vector <int> ipiv (n * nrhs, 0);

		std::vector <double> x ((ku + kl) * (ku + kl) * 4 * np * np * nrhs, 0.0);
		std::vector <int> xipiv (2 * (ku + kl) * np * nrhs, 0);

		std::vector <double> b (ldb * nrhs, 0.0);
		std::vector <double> bcopy (ldb * nrhs, 0.0);
	
		for (int i = 0; i < 10; ++i) {
			int info;

			srand (1);
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < nrhs; ++j) {
					for (int k = kl; k < 2 * kl + ku + 1; ++k) {
						matrix [j * lda * ldaa + (i + kl) * lda + k] = rand () % 100;
						matrixcopy [j * lda * ldaa + (i + kl) * lda + k] = matrix [j * lda * ldaa + (i + kl) * lda + k];
					}
					b [j * ldb + i + kl] = rand () % 100;
					bcopy [j * ldb + i + kl] = b [j * ldb + i + kl];
				}
			}
	
			if (id != 0) {
				for (int j = 0; j < nrhs; ++j) {
					for (int i = 0; i < kl; ++i) {
						for (int k = 2 * kl + ku - i; k < 2 * kl + ku + 1; ++k) {
							matrix [j * lda * ldaa + (i) * lda + k] = rand () % 100;
							matrixcopy [j * lda * ldaa + (i) * lda + k] = matrix [j * lda * ldaa + i * lda + k];
						}
					}
				}
			}
			for (int j = 0; j < nrhs; ++j) {
				for (int i = kl; i < kl + ku; ++i) {
					for (int k = kl; k < 2 * kl + ku - i; ++k) {
						matrix [j * lda * ldaa + (i) * lda + k] = 0.0;
						matrixcopy [j * lda * ldaa + (i) * lda + k] = matrix [j * lda * ldaa + i * lda + k];
					}
				}
			}

			if (id != np - 1) {
				for (int j = 0; j < nrhs; ++j) {
					for (int i = kl + n; i < kl + n + ku; ++i) {
						for (int k = kl; k < 2 * kl + ku - i + n; ++k) {
							matrix [j * lda * ldaa + (i) * lda + k] = rand () % 100;
							matrixcopy [j * lda * ldaa + (i) * lda + k] = matrix [j * lda * ldaa + i * lda + k];
						}
					}
				}
			}
			for (int j = 0; j < nrhs; ++j) {
				for (int i = n; i < kl + n; ++i) {
					for (int k = 2 * kl + ku - i + n; k < 2 * kl + ku + 1; ++k) {
						matrix [j * lda * ldaa + (i) * lda + k] = 0.0;
						matrixcopy [j * lda * ldaa + (i) * lda + k] = matrix [j * lda * ldaa + i * lda + k];
					}
				}
			}

			try {
				linalg::p_block_banded_factorize (id, np, n - ntop - nbot, kl, ku, &matrix [0], &ipiv [0], &x [0], &xipiv [0], &bufferl [0], &bufferr [0], &info, nrhs, lda, ldaa);

				linalg::p_block_banded_solve (id, np, n - ntop - nbot, kl, ku, &matrix [0], &ipiv [0], &b [kl], &x [0], &xipiv [0], &bufferl [0], &bufferr [0], &info, nrhs, lda, ldaa, ldb);
			} catch (std::exception& except) {
				std::cout << except.what () << '\n';
			}

			std::vector <double> edge_0 (nrhs * ku), edge_n (nrhs * kl), redge_0 (nrhs * kl), redge_n (nrhs * ku);

			for (int i = 0; i < nrhs; ++i) {
				for (int j = 0; j < ku; ++j) {
					edge_0 [i * ku + j] = b [i * ldb + j + kl];
				}
				for (int j = 0; j < kl; ++j) {
					edge_n [i * kl + j] = b [i * ldb + n + j];
				}
			}

			if (id != 0) {
				mess.send (nrhs * ku, &edge_0 [0], id - 1, 0);
			}
			if (id != np - 1) {
				mess.recv (nrhs * ku, &redge_n [0], id + 1, 0);
			}
			if (id != np - 1) {
				mess.send (nrhs * kl, &edge_n [0], id + 1, 1);
			}
			if (id != 0) {
				mess.recv (nrhs * kl, &redge_0 [0], id - 1, 1);
			}

			for (int i = 0; i < nrhs; ++i) {
				for (int j = 0; j < kl; ++j) {
					b [i * ldb + j] = redge_0 [i * kl + j];
				}
				for (int j = 0; j < kl; ++j) {
					b [i * ldb + n + kl + j] = redge_n [i * ku + j];
				}
			}

			for (int j = 0; j < nrhs; ++j) {
				for (int i = 0; i < n + kl + ku; ++i) {
					for (int k = -ku; k < kl + 1; ++k) {
						if (i + k >= 0 && i + k < n + kl + ku) {
							bcopy [j * ldb + k + i] -= matrixcopy [j * lda * ldaa + (i) * lda + k + kl + ku] * b [j * ldb + i];
						}
					}
				}
			}

			for (int j = 0; j < nrhs; ++j) {
				for (int i = 0; i < kl + ku + n; ++i) {
					TSM_ASSERT_DELTA ("Banded block solver failure", bcopy [j * ldb + i], 0.0, TEST_TINY);
				}
			}
		}
	}
};

