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
	io::parameters config ("../input/config.yaml");
	// io::parameters <double> params ("../input/block_parameters.txt");
	// omp_set_num_threads(params.nmp);

	int id = mess.get_id ();
	int np = mess.get_np ();
	int n = config.get <int> ("grid.z.points"), nrhs = config.get <int> ("grid.nrhs"), ntop = 0, nbot = 0;
	if (id != 0) {
		ntop = config.get <int> ("grid.nb");
	}
	if (id != np - 1) {
		nbot = config.get <int> ("grid.nb");
	}
	int lda = n + ntop + nbot, ldx = 0, ldb = lda;
	// 
	// std::vector <double> sub (n * nrhs, 0.0);
	// std::vector <double> diag (n * nrhs, 0.0);
	// std::vector <double> sup (n * nrhs, 0.0);
	// std::vector <double> subcopy (n * nrhs, 0.0);
	// std::vector <double> diagcopy (n * nrhs, 0.0);
	// std::vector <double> supcopy (n * nrhs, 0.0);
	// std::vector <double> supsup (n * nrhs, 0.0);
	// std::vector <int> ipiv (n * nrhs, 0);
	// 
	// std::vector <double> x (2 * n * nrhs + 8 * np * nrhs, 0.0);
	// std::vector <int> xipiv (2 * np * nrhs, 0);
	// 
	// std::vector <double> b (n * nrhs, 0.0);
	// std::vector <double> bcopy (n * nrhs, 0.0);
	// 
	// int info;
	// 
	// srand (1);
	// for (int i = 0; i < n; ++i) {
	// 	sub [i] = rand () % 100;
	// 	subcopy [i] = sub [i];
	// 	diag [i] = rand () % 100;
	// 	diagcopy [i] = diag [i];
	// 	sup [i] = rand () % 100;
	// 	supcopy [i] = sup [i];
	// 	b [i] = rand () % 100;
	// 	bcopy [i] = b [i];
	// 	for (int j = 1; j < nrhs; ++j) {
	// 		sub [i + j * n] = sub [i];
	// 		subcopy [i + j * n] = sub [i];
	// 		diag [i + j * n] = diag [i];
	// 		diagcopy [i + j * n] = diag [i];
	// 		sup [i + j * n] = sup [i];
	// 		supcopy [i + j * n] = sup [i];
	// 		b [i + j * n] = b [i];
	// 		bcopy [i + j * n] = b [i];
	// 	}
	// }
	// 
	// try {
	// 	utils::p_block_tridiag_factorize (id, np, n - ntop - nbot, &sub [0], &diag [0], &sup [0], &supsup [0], &ipiv [0], &x [0], &xipiv [0], &info, nrhs);
	// 	
	// 	utils::p_block_tridiag_solve (id, np, n - ntop - nbot,  &sub [0], &diag [0], &sup [0], &supsup [0], &ipiv [0], &b [0], &x [0], &xipiv [0], &info, nrhs);
	// } catch (std::exception& except) {
	// 	std::cout << except.what () << '\n';
	// }
	// 
	// for (int i = 0; i < n; ++i) {
	// 	for (int j = 0; j < nrhs; ++j) {
	// 		std::cout << "[" << id << "] " << subcopy [i] << " " << diagcopy [i] << " " << supcopy [i] << " = " << bcopy [i] << " => " << sub [i] << " " << diag [i] << " " << sup [i] << " = " << b [i] << " | ";
	// 	}
	// 	std::cout << "\n";
	// }
	// 
	// for (int j = 0; j < nrhs; ++j) {
	// 	bcopy [j * n] -= diagcopy [j * n] * b [j * n] + supcopy [j * n] * b [1 + j * n];
	// 	for (int i = 1; i < n - 1; ++i) {
	// 		bcopy [i + j * n] -= subcopy [i + j * n] * b [i - 1 + j * n] + diagcopy [i + j * n] * b [i + j * n] + supcopy [i + j * n] * b [i + 1 + j * n];
	// 	}
	// 	bcopy [n - 1 + j * n] -= subcopy [n - 1 + j * n] * b [n - 2 + j * n] + diagcopy [n - 1 + j * n] * b [n - 1 + j * n];
	// }
	// 
	// double above;
	// double below;
	// 
	// if (id != 0) {
	// 	mess.send (1, &b [0], id - 1, 0);
	// }
	// if (id != np - 1) {
	// 	mess.recv (1, &below, id + 1, 0);
	// 	bcopy [n - 1] -= supcopy [n - 1] * below;
	// }
	// if (id != np - 1) {
	// 	mess.send (1, &b [n - 1], id + 1, 1);
	// }
	// if (id != 0) {
	// 	mess.recv (1, &above, id - 1, 1);
	// 	bcopy [0] -= subcopy [0] * above;
	// }
	// 
	// for (int i = 0; i < n; ++i) {
	// 	printf ("FINAL!: [%i] ", id);
	// 	for (int j = 0; j < nrhs; ++j) {
	// 		printf ("%f ", bcopy [i + j * n]);
	// 	}
	// 	printf ("\n");
	// }
	
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

