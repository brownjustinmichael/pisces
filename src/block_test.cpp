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
#include "utils/block_solver.hpp"
#include "utils/utils.hpp"
#include <omp.h>

int main (int argc, char *argv[])
{
	bases::messenger mess (&argc, &argv, 2);
	io::parameters <double> params ("../input/block_parameters.txt");
	omp_set_num_threads(params.nmp);

	int n = params.n, nrhs = params.nrhs, ntimes = params.timesteps, ntop = 0, nbot = 0;
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
	
	double begin = omp_get_wtime ();

	utils::p_block_matrix_factorize (&mess, n, ntop, nbot, &a [0], &ipiv [0], &x [0], &xipiv [0], &ns [0], &info, lda, ldx);

	double mid = omp_get_wtime ();

	for (int i = 0; i < ntimes; ++i) {
		utils::p_block_matrix_solve (&mess, n, ntop, nbot, &a [0], &ipiv [0], &b [0], &x [0], &xipiv [0], &ns [0], &info, nrhs, lda, ldx, ldb);
	}
	
	double end = omp_get_wtime ();
	
	// utils::matrix_matrix_multiply (lda, nrhs, lda, -1.0, &acopy [0], &b [0], 1.0, &bcopy [0]);
	
	// for (int i = 0; i < lda; ++i) {
	// 	printf ("%d ", mess.get_id ());
	// 	for (int j = 0; j < lda; ++j) {
	// 		printf ("%f ", a [j * lda + i]);
	// 	}
	// 	printf ("= ");
	// 	for (int j = 0; j < nrhs; ++j) {
	// 		printf ("%f ", b [j * ldb + i]);
	// 	}
	// 	printf ("\n");
	// }
	
	// for (int i = 0; i < lda; ++i) {
	// 	printf ("r [%d]:", mess.get_id ());
	// 	for (int j = 0; j < nrhs; ++j) {
	// 		if (!(bcopy [j * ldb + i] < 1.0e-6 && bcopy [j * ldb + i] > -1.0e-6)) {
	// 			printf ("%f ", bcopy [j * ldb + i]);
	// 		}
	// 	}
	// 	printf ("\n");
	// }
	
	printf ("[%d]: Time: %f + %f = %f\n", mess.get_id (), (mid - begin), (end - mid), (end - begin));

	// int nm = 5, nbtot = 0, info, ntot = 0, ncur = 0;
	// std::vector <int> n (nm);
	// std::vector <int> nb (nm);
	// 
	// for (int i = 0; i < nm; ++i) {
	// 	n [i] = 8;
	// 	if (i != nm - 1) {
	// 		nb [i] = i > 2 ? 1 : 3;
	// 	} else {
	// 		nb [i] = 0;
	// 	}
	// 	nbtot += nb [i];
	// 	ntot += n [i] - nb [i];
	// }
	// 
	// std::vector <double *> a (nm);
	// std::vector <double *> b (nm);
	// std::vector <int *> ipiv (nm);
	// 
	// std::vector <double> x (nbtot * nbtot);
	// std::vector <int> xipiv (nbtot);
	// 
	// std::vector <std::vector <double>> ta (nm);
	// std::vector <std::vector <double>> tb (nm);
	// std::vector <std::vector <int>> tipiv (nm);
	// 
	// for (int i = 0; i < nm; ++i) {
	// 	ta [i].resize (n [i] * n [i]);
	// 	tb [i].resize (n [i]);
	// 	tipiv [i].resize (n [i]);
	// 	
	// 	a [i] = &(ta [i] [0]);
	// 	b [i] = &(tb [i] [0]);
	// 	ipiv [i] = &(tipiv [i] [0]);
	// }
	// 
	// std::vector <double> matrix (ntot * ntot, 0.0);
	// std::vector <double> rhs (ntot, 0.0);
	// std::vector <double> lhs (ntot, 0.0);
	// std::vector <int> mipiv (ntot, 0);
	// 
	// srand (2);
	// for (int i = 0; i < nm; ++i) {
	// 	for (int j = 0; j < n [i]; ++j) {
	// 		tb [i] [j] = rand () % 100;
	// 		printf ("%f = ", b [i] [j]);
	// 		rhs [j + ncur] += tb [i] [j];
	// 		for (int k = 0; k < n [i]; ++k) {
	// 			ta [i] [k * n [i] + j] = rand () %100;
	// 			printf ("%f, ", a [i] [k * n [i] + j]);
	// 			matrix [(k + ncur) * ntot + (j + ncur)] += ta [i] [k * n [i] + j];
	// 		}
	// 		printf ("\n");
	// 	}
	// 	ncur += n[i] - nb [i];
	// }
	// 
	// for (int j = 0; j < ntot; ++j) {
	// 	printf ("%f = ", rhs [j]);
	// 	for (int k = 0; k < ntot; ++k) {
	// 		printf ("%f, ", matrix [k * ntot + j]);
	// 	}
	// 	printf ("\n");
	// }
	// // 
	// // utils::matrix_factorize (ntot, ntot, &matrix [0], &mipiv [0]);
	// // utils::matrix_solve (ntot, &matrix [0], &mipiv [0], &rhs [0]);
	// 
	// printf ("Factorizing...\n");
	// utils::block_matrix_factorize (nm, &n [0], &nb [0], &a [0], &ipiv [0], &x [0], &xipiv [0], &info);
	// 
	// for (int i = 0; i < nm; ++i) {
	// 	for (int j = 0; j < n [i]; ++j) {
	// 		printf ("%f = ", b [i] [j]);
	// 		for (int k = 0; k < n [i]; ++k) {
	// 			printf ("%f, ", a [i] [k * n [i] + j]);
	// 		}
	// 		printf ("\n");
	// 	}
	// }
	// 
	// printf ("Solving...\n");
	// utils::block_matrix_solve (nm, &n [0], &nb [0], &a [0], &ipiv [0], &x [0], &xipiv [0], &b [0], &info);
	// printf ("Done.\n");
	// 
	// ncur = 0;
	// for (int i = 0; i < nm; ++i) {
	// 	for (int j = 0; j < n [i] - nb [i]; ++j) {
	// 		lhs [ncur] = b [i] [j];
	// 		++ncur;
	// 	}
	// }
	// 
	// for (int i = 0; i < nm; ++i) {
	// 	for (int j = 0; j < n [i]; ++j) {
	// 		printf ("B [%d]: %f\n", i, b [i] [j]);
	// 	}
	// }
	// 
	// for (int i = 0; i < ntot; ++i) {
	// 	printf ("LHS: %f\n", lhs [i]);
	// }
	// 
	// utils::matrix_matrix_multiply (ntot, 1, ntot, 1.0, &matrix [0], &lhs [0], -1.0, &rhs [0]);
	// 
	// for (int i = 0; i < ntot; ++i) {
	// 	printf ("%f\n", rhs [i]);
	// }
	
	return 0;
}

