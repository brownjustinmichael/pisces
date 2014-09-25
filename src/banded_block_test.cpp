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
	std::stringstream debug;
	std::stringstream matrix_debug;
	
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
	
	
	if (argc <= 1) {
		config_filename = "../input/config.yaml";
	} else {
		config_filename = argv [1];
	}
	io::parameters config (config_filename);

	// io::parameters <double> params ("../input/block_parameters.txt");
	// omp_set_num_threads(params.nmp);

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

		for (int j = 0; j < nrhs; ++j) {
			for (int i = -kl; i < n + ku; ++i) {
				for (int k = 0; k < lda; ++k) {
					debug << matrix [j * lda * ldaa + (i + kl) * lda + k] << " ";
				}
				DEBUG ("[" << id << "] " << debug.str () << "= " << ((i >= 0 && i < n) ? b [j * ldb + i + kl] : 0));
				debug.str ("");
			}
		}
	
		try {
			utils::p_block_banded_factorize (id, np, n - ntop - nbot, kl, ku, &matrix [0], &ipiv [0], &x [0], &xipiv [0], &bufferl [0], &bufferr [0], &info, nrhs, lda, ldaa);

			utils::p_block_banded_solve (id, np, n - ntop - nbot, kl, ku, &matrix [0], &ipiv [0], &b [kl], &x [0], &xipiv [0], &bufferl [0], &bufferr [0], &info, nrhs, lda, ldaa, ldb);
		} catch (std::exception& except) {
			std::cout << except.what () << '\n';
		}

		for (int j = 0; j < nrhs; ++j) {
			for (int i = 0; i < n + ku + kl; ++i) {
				for (int k = 0; k < lda; ++k) {
					debug << matrix [j * lda * ldaa + i * lda + k] << " ";
				}
				DEBUG ("[" << id << "] " << debug.str ());
				debug.str ("");
			}
		}

		for (int j = 0; j < nrhs; ++j) {
			DEBUG ("Mid " << j);
			for (int i = 0; i < n; ++i) {
				for (int k = 0; k < 2 * kl + ku + 1; ++k) {
					matrix_debug << matrixcopy [j * lda * ldaa + (i + kl) * lda + k] << " ";
				}
				debug << bcopy [j * ldb + i + kl] << " ";
				DEBUG (matrix_debug.str () << " = " << debug.str ());
				matrix_debug.str ("");
				debug.str ("");
			}
		}

		std::vector <double> edge_0 (nrhs * ku), edge_n (nrhs * kl), redge_0 (nrhs * kl), redge_n (nrhs * ku);

		for (int j = 0; j < nrhs; ++j) {
			DEBUG ("BEFORE SEND " << j);
			for (int i = 0; i < kl + ku + n; ++i) {
				for (int k = 0; k < 2 * kl + ku + 1; ++k) {
					matrix_debug << matrixcopy [j * lda * ldaa + (i) * lda + k] << " ";
				}
				debug << b [j * ldb + i] << " ";
				DEBUG (matrix_debug.str () << " = " << debug.str ());
				matrix_debug.str ("");
				debug.str ("");
			}
		}

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
			DEBUG ("AFTER RECV " << j);
			for (int i = 0; i < kl + ku + n; ++i) {
				for (int k = 0; k < 2 * kl + ku + 1; ++k) {
					matrix_debug << matrixcopy [j * lda * ldaa + (i) * lda + k] << " ";
				}
				debug << b [j * ldb + i] << " ";
				DEBUG (matrix_debug.str () << " = " << debug.str () << "=> " << bcopy [j * ldb + i]);
				matrix_debug.str ("");
				debug.str ("");
			}
		}

		for (int j = 0; j < nrhs; ++j) {
			for (int i = 0; i < n + kl + ku; ++i) {
				for (int k = -ku; k < kl + 1; ++k) {
					if (i + k >= 0 && i + k < n + kl + ku) {
						DEBUG (k << " " << i << " " << j * ldb + k + i << " " << bcopy [j * ldb + i + k] << " " << matrixcopy [j * lda * ldaa + (i) * lda + k + kl + ku] << " " << b [j * ldb + i]);
						bcopy [j * ldb + k + i] -= matrixcopy [j * lda * ldaa + (i) * lda + k + kl + ku] * b [j * ldb + i];
					}
				}
			}
		}

		for (int j = 0; j < nrhs; ++j) {
			DEBUG ("END " << j);
			for (int i = 0; i < kl + ku + n; ++i) {
				for (int k = 0; k < 2 * kl + ku + 1; ++k) {
					matrix_debug << matrix [j * lda * ldaa + (i) * lda + k] << " ";
				}
				debug << bcopy [j * ldb + i] << " ";
				if (bcopy [j * ldb + i] > 1.0e-8) {
					ERROR ("FAILED TO ACHIEVE TOLERANCE");
				}
				DEBUG (matrix_debug.str () << " = " << debug.str ());
				matrix_debug.str ("");
				debug.str ("");
			}
		}
	}

	return 0;
}

