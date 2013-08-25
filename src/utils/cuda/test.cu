/*!**********************************************************************
 * \file main_cuda.cpp
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-14.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

// #include "config.hpp"
// #include "one_d/cuda/element_one_d_cuda.hpp"
#include <vector>
#include <stdio.h>
#include "dgetrs.h"
#include <stdlib.h>

extern "C" {void dgetrf_ (int* n, int* m, double* a, int* lda, int* ipiv, int* info);}

int main (int argc, char *argv[])
{	
	// int id = 0, n_elements = 1;
	// 
	// // Initialize messenger
	// bases::messenger <double> process_messenger (&argc, &argv, 2);
	// 
	// id = process_messenger.get_id ();
	// n_elements = process_messenger.get_np ();
	// 
	// log_config::update_name (id);
	// 
	// // The program runs through the execution flags.
	// while ((argc > 1) && (argv [1] [0] == '-')) {
	// 	switch (argv [1] [1]) {
	// 		// Debug switch
	// 		case 'D':
	// 			log_config::update_severity (atoi (&(argv [1] [2])));
	// 			break;
	// 	}
	// 	--argc;
	// 	++argv;
	// }
	// 
	// TRACE ("Command line arguments read, beginning setup.");
	// 
	// 	
	// io::parameter_map inputParams;
	// io::read_params_txt parameters ("../input/parameters.txt");
	// inputParams = parameters.load_params();
	// 
	// 
	// int n = inputParams ["gridpoints"].asInt / n_elements;
	// double position_0 = -1.0 + 2.0 / n_elements * id;
	// double position_n = -1.0 + 2.0 / n_elements * (id + 1);
	// int excess_0;
	// int excess_n;
	// if (id == 0) {
	// 	excess_0 = 0;
	// } else {
	// 	excess_0 = 1;
	// }
	// if (id == n_elements - 1) {
	// 	excess_n = 0;
	// } else {
	// 	excess_n = 1;
	// }
	// int name = id;
	// 
	// if (id != 0) {
	// 	TRACE ("Adding boundary to " << name << " at 0 at processor " << id - 1);
	// 	process_messenger.add_boundary (one_d::edge_0, id - 1);
	// }
	// if (id != n_elements - 1) {
	// 	TRACE ("Adding boundary to " << name << " at n - 1 at processor " << id + 1);
	// 	process_messenger.add_boundary (one_d::edge_n, id + 1);
	// }
	// 
	// one_d::chebyshev::cuda::fft_element <double> element (n, position_0, position_n, excess_0, excess_n, name, inputParams, &process_messenger, 0x00);
	// 
	// element.setup ();
	// 
	// try {
	// 	element.run ();
	// } catch (...) {
	// 	FATAL ("Fatal error occurred. Check log.");
	// 	return 1;
	// }
	
	// INFO ("Main complete.");

	int n = 1024;
	int nrhs = 6;
	int lda = n + 106, ldb = n + 20;
	int info;
	double range = 10.0;
	
	double* a_dev, *b_dev;
	int* ipiv_dev;
	std::vector <double> matrix (n * lda);
	std::vector <double> mcopy (n * lda);
	std::vector <double> x (n * ldb), xcopy (n * ldb);
	std::vector <int> ipiv (n);
	
	for (int i = 0; i < n; ++i) {
		// printf ("Before: ");
		for (int j = 0; j < n; ++j) {
			matrix [i + j * lda] = (double) rand () / (double) RAND_MAX * 2.0 * range - range;
			mcopy [i + j * lda] = matrix [i + j * lda];
			// printf ("%e ", matrix [i + j * lda]);
		}
		// printf ("= ");
		for (int j = 0; j < nrhs; ++j) {
			x [i + j * ldb] = (double) rand () / (double) RAND_MAX * 2.0 * range - range;
			xcopy [i + j * ldb] = x [i + j * ldb];
			// printf ("%e ", x[i + j * ldb]);
		}
		// printf ("\n");
	}
	
	
	dgetrf_ (&n, &n, &matrix [0], &lda, &ipiv [0], &info);
	
	// solve_lu (n, nrhs, &matrix [0], lda, &ipiv [0], &x [0], ldb, &info);
		
	cudaMalloc (&a_dev, sizeof (double) * n * lda);
	cudaMalloc (&b_dev, sizeof (double) * nrhs * ldb);
	cudaMalloc (&ipiv_dev, sizeof (int) * n);
	
	cudaMemcpy (a_dev, &matrix [0], sizeof (double) * n * lda, cudaMemcpyHostToDevice);
	cudaMemcpy (b_dev, &x [0], sizeof (double) * nrhs * ldb, cudaMemcpyHostToDevice);
	cudaMemcpy (ipiv_dev, &ipiv [0], sizeof (int) * n, cudaMemcpyHostToDevice);
	
	solve_lu (n, nrhs, a_dev, lda, ipiv_dev, b_dev, ldb, &info);
	
	cudaMemcpy (&x [0], b_dev, sizeof (double) * nrhs * ldb, cudaMemcpyDeviceToHost);
	
	cudaDeviceSynchronize ();
	
	cudaFree (ipiv_dev);
	cudaFree (a_dev);
	cudaFree (b_dev);
	
	if (info != 0) {
		printf ("Error");
	}
	
	for (int i = 0; i < n; ++i) {
	// 	printf ("After: ");
	// 	for (int j = 0; j < n; ++j) {
	// 		printf ("%f ", matrix [i + j * n]);
	// 	}
		// printf ("= ");
		// for (int j = 0; j < nrhs; ++j) {
		// 	printf ("%f ", x[i + j * ldb]);
		// }
		// printf ("\n");
	}
	
	for (int i = 0; i < n; ++i) {
		for (int k = 0; k < nrhs; ++k) {
			for (int j = 0; j < n; ++j) {
				xcopy [i + k * ldb] -= mcopy [i + j * lda] * x [j + k * ldb];
			}
			if (abs (xcopy [i + k * ldb]) > 1.e-4) {
				printf ("Problem %e.\n", xcopy [i + k * ldb]);
			}
		}
	}
	
	printf ("Done\n");
	
	return 0;
}