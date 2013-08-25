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
#include "utils/solver_utils.hpp"
#include "utils/cuda/dgetrs.h"

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
	// one_d::chebyshev::cuda::fft_element <double> element (n, excess_0, position_0, excess_n, position_n, name, inputParams, &process_messenger, 0x00);
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

	int n = 6;
	int info;
	std::vector <double> matrix (n * n);
	std::vector <double> x (n);
	std::vector <int> ipiv (n);
	
	for (int i = 0; i < n; ++i) {
		printf ("Before: ");
		for (int j = 0; j < n; ++j) {
			matrix [i + j * n] = i + j;
			printf ("%f ", matrix [i + j * n]);
		}
		x [i] = i;
		printf ("= %f\n", x[i]);
	}
	
	char charN = 'N';
	int ione = 1;

	// utils::matrix_factorize (n, n, &matrix [0], &ipiv [0], &info);
	// dgetrscuda_ (&charN, &n, &ione, &matrix [0], &n, &ipiv [0], &x [0], &n, &info);
	
	if (info != 0) {
		printf ("Error");
	}
	
	for (int i = 0; i < n; ++i) {
		printf ("After: ");
		for (int j = 0; j < n; ++j) {
			printf ("%f ", matrix [i + j * n]);
		}
		printf ("= %f\n", x[i]);
	}
	
	return 0;
}