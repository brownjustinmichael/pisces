/*!**********************************************************************
 * \file block_test.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-11-05.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cxxtest/TestSuite.h>

#include "messenger/messenger.hpp"
#include "io/io.hpp"
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include "plan-solver/block_solver.hpp"
#include "linalg/linalg.hpp"
#include "linalg/utils.hpp"
#include <iostream>

#define TEST_TINY 1.e-4

class mpi_solver_test_suite : public CxxTest::TestSuite
{
public:
	void test_block_solver () {
		utils::messenger mess;
		
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
	
		utils::matrix_copy (lda, lda, &a [0], &acopy [0], lda, lda);
	
		clock_t corig_begin, corig_end;
		std::chrono::time_point <std::chrono::system_clock> orig_begin, orig_end;
	
		corig_begin = clock ();
		orig_begin = std::chrono::system_clock::now ();
	
		utils::matrix_factorize (n, n, &a [0], &ipiv [0], &info, lda);
		utils::matrix_solve (n, &a [0], &ipiv [0], &b [0], &info, nrhs, lda, ldb);
	
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
		
			utils::p_block_matrix_factorize (mess.get_id (), mess.get_np (), n, ntop, nbot, &a [0], &ipiv [0], &x [0], &xipiv [0], &ns [0], &info, lda, ldx);
	
			cmid = clock ();
			mid = std::chrono::system_clock::now ();
		
			mb += mid - begin;
			cmb += cmid - cbegin;
		
			utils::p_block_matrix_solve (mess.get_id (), mess.get_np (), n, ntop, nbot, &a [0], &ipiv [0], &b [0], &x [0], &xipiv [0], &ns [0], &info, nrhs, lda, ldx, ldb);
	
			cend = clock ();
			end = std::chrono::system_clock::now ();
		
			em += end - mid;
			eb += end - begin;
			cem += cend - cmid;
			ceb += cend - cbegin;
			
			utils::matrix_matrix_multiply (n + ntop + nbot, nrhs, n + ntop + nbot, 1.0, &acopy [0], &b [0], -1.0, &bcopy [0]);
			
			if (id != 0) {
				utils::matrix_copy (ntop, nrhs, &bcopy [0], &bufftop [0], ldb, ntop);
				mess.send (ntop * nrhs, &bufftop [0], id - 1, 0);
			}
			if (id != np - 1) {
				utils::matrix_copy (nbot, nrhs, &bcopy [n + ntop], &buffbot [0], ldb, nbot);
				mess.recv (nbot * nrhs, &rbuffbot [0], id + 1, 0);
				mess.send (nbot * nrhs, &buffbot [0], id + 1, 1);
				utils::matrix_add_scaled (nbot, nrhs, 1.0, &rbuffbot [0], &bcopy [n + ntop], nbot, ldb);
			}
			if (id != 0) {
				mess.recv (ntop * nrhs, &rbufftop [0], id - 1, 1);
				utils::matrix_add_scaled (ntop, nrhs, 1.0, &rbufftop [0], &bcopy [0], ntop, ldb);
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
};

