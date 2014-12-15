/*!**********************************************************************
 * \file block_test.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-11-05.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cxxtest/TestSuite.h>

#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include "linalg/linalg.hpp"
#include <omp.h>

#define TEST_TINY 1.e-6

class linalg_test_suite : public CxxTest::TestSuite
{
public:
	void test_matrix_solver () {
		int n = 1000, nrhs = 10000;
		int lda = n, ldx = 0, ldb = lda;
		
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
	
		// omp_set_num_threads (8);
	
		corig_begin = clock ();
		orig_begin = std::chrono::system_clock::now ();
	
		linalg::matrix_factorize (n, n, &a [0], &ipiv [0], &info, lda);
		linalg::matrix_solve (n, &a [0], &ipiv [0], &b [0], &info, nrhs, lda, ldb);
	
		corig_end = clock ();
		orig_end = std::chrono::system_clock::now ();
		
	
		std::chrono::duration <double> orig_total = orig_end - orig_begin;
	
		printf ("Orig CPU Time: %f\n", ((double) (corig_end - corig_begin))/CLOCKS_PER_SEC);
		printf ("Orig Wall Time: %f\n", (double) orig_total.count ());
	}
};

