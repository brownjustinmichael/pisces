/*!**********************************************************************
 * \file block_test.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-11-05.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "utils/messenger.hpp"
#include "io/io.hpp"
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include "utils/block_solver.hpp"
#include "utils/solver_utils.hpp"
#include "utils/utils.hpp"
#include <sstream>

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
		for (int j = 0; j < nrhs; ++j) {
			sub [i * nrhs + j] = rand () % 100;
			subcopy [i * nrhs + j] = sub [i * nrhs + j];
			diag [i * nrhs + j] = rand () % 100;
			diagcopy [i * nrhs + j] = diag [i * nrhs + j];
			sup [i * nrhs + j] = rand () % 100;
			supcopy [i * nrhs + j] = sup [i * nrhs + j];
			b [i * nrhs + j] = rand () % 100;
			bcopy [i * nrhs + j] = b [i * nrhs + j];
		}
	}
	
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < nrhs; ++j) {
			DEBUG ("[" << id << "] " << subcopy [i] << " " << diagcopy [i] << " " << supcopy [i] << " = " << bcopy [i]);
		}
	}

	try {
		utils::p_block_tridiag_factorize (id, np, n - ntop - nbot, &sub [0], &diag [0], &sup [0], &supsup [0], &ipiv [0], &x [0], &xipiv [0], &info, nrhs);
		// utils::p_block_tridiag_factorize (id, np, n - ntop - nbot, &sub [0], &diag [0], &sup [0], &supsup [0], &ipiv [0], &x [0], &xipiv [0], &info, 1);

		utils::p_block_tridiag_solve (id, np, n - ntop - nbot,  &sub [0], &diag [0], &sup [0], &supsup [0], &ipiv [0], &b [0], &x [0], &xipiv [0], &info, nrhs);
		// utils::p_block_tridiag_solve (id, np, n - ntop - nbot,  &sub [0], &diag [0], &sup [0], &supsup [0], &ipiv [0], &b [0], &x [0], &xipiv [0], &info, 1, -1, -1, nrhs);
		// utils::p_block_direct_tridiag_solve (id, np, n - ntop - nbot,  &sub [0], &diag [0], &sup [0], &supsup [0], &b [0], nrhs, n);
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
	
	std::vector <double> edge_0 (nrhs), edge_n (nrhs), redge_0 (nrhs), redge_n (nrhs);
	
	for (int i = 0; i < nrhs; ++i) {
		edge_0 [i] = b [i * nrhs];
		edge_n [i] = b [i * nrhs + n - 1];
	}

	if (id != 0) {
		// mess.send (1, &b [0], id - 1, 0);
		mess.send (nrhs, &edge_0 [0], id - 1, 0);
	}
	if (id != np - 1) {
		mess.recv (nrhs, &redge_n [0], id + 1, 0);
		// mess.recv (1, &below, id + 1, 0);
		// bcopy [n - 1] -= supcopy [n - 1] * below;
		for (int i = 0; i < nrhs; ++i) {
			bcopy [i * nrhs + n - 1] -= supcopy [i * nrhs + n - 1] * redge_n [i];
		}
	}
	if (id != np - 1) {
		// mess.send (1, &b [n - 1], id + 1, 1);
		mess.send (nrhs, &edge_n [0], id + 1, 1);
	}
	if (id != 0) {
		mess.recv (nrhs, &redge_0 [0], id - 1, 1);
		// mess.recv (1, &above, id - 1, 1);
		// bcopy [0] -= subcopy [0] * above;
		for (int i = 0; i < nrhs; ++i) {
			bcopy [i * nrhs] -= subcopy [i * nrhs] * redge_0 [i];
		}
	}

	std::stringstream sub_debug, diag_debug, sup_debug, debug;
	for (int j = 0; j < nrhs; ++j) {
		for (int i = 0; i < n; ++i) {
			sub_debug << subcopy [i * nrhs + j] << " ";
			diag_debug << diagcopy [i * nrhs + j] << " ";
			sup_debug << supcopy [i * nrhs + j] << " ";
			debug << bcopy [i * nrhs + j] << " ";
		}
		DEBUG (sub_debug.str () << " | " << diag_debug.str () << " | " << sup_debug.str () << " = " << debug.str ());
		sub_debug.str ("");
		diag_debug.str ("");
		sup_debug.str ("");
		debug.str ("");
	}
	return 0;
}

