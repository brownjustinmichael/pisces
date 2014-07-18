/*!**********************************************************************
 * \file block_test.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-11-05.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "utils/messenger.hpp"
#include "utils/io.hpp"
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include "utils/block_solver.hpp"
#include "utils/solver_utils.hpp"
#include "utils/utils.hpp"

int main (int argc, char *argv[])
{
	utils::messenger mess (&argc, &argv);
	log_config::configure (&argc, &argv, mess.get_id (), "process_%d.log");
	std::string config_filename;
	
	if (argc <= 1) {
		config_filename = "../input/config.yaml";
	} else {
		config_filename = argv [1];
	}
	io::parameters config (config_filename);

	// io::parameters <double> params ("../input/block_parameters.txt");
	// omp_set_num_threads(params.nmp);

	int id = mess.get_id ();
	int np = mess.get_np ();
	int n = config.get <int> ("grid.z.points"), nrhs = config.get <int> ("grid.x.points"), ntop = 0, nbot = 0;
	if (id != 0) {
		ntop = 1;
	}
	if (id != np - 1) {
		nbot = 1;
	}
	int lda = n + ntop + nbot, ldx = 0, ldb = lda;

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
	for (int i = 0; i < n; ++i) {
		sub [i] = rand () % 100;
		subcopy [i] = sub [i];
		diag [i] = rand () % 100;
		diagcopy [i] = diag [i];
		sup [i] = rand () % 100;
		supcopy [i] = sup [i];
		b [i] = rand () % 100;
		bcopy [i] = b [i];
		for (int j = 1; j < nrhs; ++j) {
			sub [i + j * n] = sub [i];
			subcopy [i + j * n] = sub [i];
			diag [i + j * n] = diag [i];
			diagcopy [i + j * n] = diag [i];
			sup [i + j * n] = sup [i];
			supcopy [i + j * n] = sup [i];
			b [i + j * n] = b [i];
			bcopy [i + j * n] = b [i];
		}
	}
	
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < nrhs; ++j) {
			DEBUG ("[" << id << "] " << subcopy [i] << " " << diagcopy [i] << " " << supcopy [i] << " = " << bcopy [i]);
		}
	}

	try {
		// utils::p_block_tridiag_factorize (id, np, n - ntop - nbot, &sub [0], &diag [0], &sup [0], &supsup [0], &ipiv [0], &x [0], &xipiv [0], &info, nrhs);
		// utils::p_block_tridiag_factorize (id, np, n - ntop - nbot, &sub [0], &diag [0], &sup [0], &supsup [0], &ipiv [0], &x [0], &xipiv [0], &info, 1);

		// utils::p_block_tridiag_solve (id, np, n - ntop - nbot,  &sub [0], &diag [0], &sup [0], &supsup [0], &ipiv [0], &b [0], &x [0], &xipiv [0], &info, nrhs);
		// utils::p_block_tridiag_solve (id, np, n - ntop - nbot,  &sub [0], &diag [0], &sup [0], &supsup [0], &ipiv [0], &b [0], &x [0], &xipiv [0], &info, 1, -1, -1, nrhs);
		utils::p_block_direct_tridiag_solve (id, np, n - ntop - nbot,  &sub [0], &diag [0], &sup [0], &supsup [0], &b [0], nrhs, n);
	} catch (std::exception& except) {
		std::cout << except.what () << '\n';
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < nrhs; ++j) {
			DEBUG ("[" << id << "] " << subcopy [i] << " " << diagcopy [i] << " " << supcopy [i] << " = " << bcopy [i] << " => " << sub [i] << " " << diag [i] << " " << sup [i] << " = " << b [i] << " | ");
		}
	}

	for (int j = 0; j < nrhs; ++j) {
		bcopy [j * n] -= diagcopy [j * n] * b [j * n] + supcopy [j * n] * b [1 + j * n];
		for (int i = 1; i < n - 1; ++i) {
			bcopy [i + j * n] -= subcopy [i + j * n] * b [i - 1 + j * n] + diagcopy [i + j * n] * b [i + j * n] + supcopy [i + j * n] * b [i + 1 + j * n];
		}
		bcopy [n - 1 + j * n] -= subcopy [n - 1 + j * n] * b [n - 2 + j * n] + diagcopy [n - 1 + j * n] * b [n - 1 + j * n];
	}

	double above;
	double below;

	if (id != 0) {
		mess.send (1, &b [0], id - 1, 0);
	}
	if (id != np - 1) {
		mess.recv (1, &below, id + 1, 0);
		bcopy [n - 1] -= supcopy [n - 1] * below;
	}
	if (id != np - 1) {
		mess.send (1, &b [n - 1], id + 1, 1);
	}
	if (id != 0) {
		mess.recv (1, &above, id - 1, 1);
		bcopy [0] -= subcopy [0] * above;
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < nrhs; ++j) {
			DEBUG ("[" << id << "] " << subcopy [i] << " " << diagcopy [i] << " " << supcopy [i] << " = " << bcopy [i]);
		}
	}
	return 0;
}

