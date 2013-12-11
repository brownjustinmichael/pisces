/*!**********************************************************************
 * \file block_test.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-11-05.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "bases/messenger.hpp"
#include "utils/io.hpp"
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include "utils/block_solver.hpp"
#include "utils/utils.hpp"

int main (int argc, char *argv[])
{
	bases::messenger mess (&argc, &argv, 2);
	io::parameters <double> params ("../input/block_parameters.txt");
	// omp_set_num_threads(params.nmp);

	int n = params.n, nrhs = params.nrhs, ntop = 0, nbot = 0;
	if (mess.get_id () != 0) {
		ntop = params.nb;
	}
	if (mess.get_id () != mess.get_np () - 1) {
		nbot = params.nb;
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
	
	std::vector <double> a (lda * lda), acopy (lda * lda);
	std::vector <double> b (lda * nrhs), bcopy (lda * nrhs);
	std::vector <int> ipiv (lda);
	std::vector <double> x (ldx * ldx);
	std::vector <int> xipiv (ldx);
	int info;
	
	srand (2);
	
	for (int i = 0; i < lda; ++i) {
		// printf ("%d ", mess.get_id ());
		for (int j = 0; j < lda; ++j) {
			a [j * lda + i] = rand () % 100;
			// printf ("%f ", a [j * lda + i]);
		}
		// printf ("= ");
		for (int j = 0; j < nrhs; ++j) {
			b [j * ldb + i] = rand () % 100;
			bcopy [j * ldb + i] = b [j * ldb + i];
			// printf ("%f ", b [j * ldb + i]);
		}
		// printf ("\n");
	}
	
	utils::matrix_copy (lda, lda, &a [0], &acopy [0], lda, lda);
	
	clock_t cbegin, cmid, cend;
	std::chrono::time_point <std::chrono::system_clock> begin, mid, end;
	
	cbegin = clock ();
	begin = std::chrono::system_clock::now ();

	utils::p_block_matrix_factorize (mess.get_id (), mess.get_np (), n, ntop, nbot, &a [0], &ipiv [0], &x [0], &xipiv [0], &ns [0], &info, lda, ldx);

	cmid = clock ();
	mid = std::chrono::system_clock::now ();
	
	utils::p_block_matrix_solve (mess.get_id (), mess.get_np (), n, ntop, nbot, &a [0], &ipiv [0], &b [0], &x [0], &xipiv [0], &ns [0], &info, nrhs, lda, ldx, ldb);
	
	cend = clock ();
	end = std::chrono::system_clock::now ();
	
	utils::matrix_matrix_multiply (n + ntop + nbot, nrhs, n + ntop + nbot, 1.0, &acopy [0], &b [0], -1.0, &bcopy [0]);
	
	// for (int i = 0; i < lda; ++i) {
	// 	printf ("%d ", mess.get_id ());
	// 	for (int j = 0; j < nrhs; ++j) {
	// 		printf ("%f ", bcopy [j * ldb + i]);
	// 	}
	// 	printf ("\n");
	// }
	
	std::chrono::duration <double> mb = mid - begin;
	std::chrono::duration <double> em = end - mid;
	std::chrono::duration <double> eb = end - begin;
	printf ("CPU Time: %f + %f = %f\n", ((double) (cmid - cbegin))/CLOCKS_PER_SEC, ((double) (cend - cmid))/CLOCKS_PER_SEC, ((double) (cend - cbegin))/CLOCKS_PER_SEC);
	printf ("Wall Time: %f + %f = %f\n", (double) mb.count (), (double) em.count (), (double) eb.count ());
	return 0;
}

