/*!**********************************************************************
 * \file block_solver.cpp
 * /Users/justinbrown/Dropbox/pisces/src/utils
 * 
 * Created by Justin Brown on 2013-11-17.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <vector>
#include "../bases/messenger.hpp"
#include "../utils/utils.hpp"
#include "../utils/solver_utils.hpp"
#include "block_solver.hpp"

namespace utils
{
	void p_block_matrix_factorize (bases::messenger* messenger_ptr, int n, int ntop, int nbot, double* a, int* ipiv, double *x, int *xipiv, int* ns, int *info, int lda, int ldx) {
		int ntot = 0, n2tot = 0, ldb = ntop + nbot, id = messenger_ptr->get_id (), np = messenger_ptr->get_np ();
		double* am = a + ntop * (lda + 1);
		double* at = a + ntop * lda;
		double* ab = a + ntop * lda + ntop + n;
		double* ar = a + (ntop + n) * lda + ntop;
		double* al = a + ntop;
		if (lda == -1) {
			lda = n + ntop + nbot;
		}
		
		printf ("Beginning...");
				
		std::vector <int> ns2;
		if (id == 0) {
			ns2.resize (np);
			for (int i = 0; i < messenger_ptr->get_np () - 1; ++i) {
				ntot += ns [i];
				ns2 [i] = (ns [i] + ns [i + 1]) * (ns [i] + ns [i + 1]);
				n2tot += ns2 [i];
			}
			ntot += ns [np - 1];
			ns2 [np - 1] = ns [np - 1] * ns [np - 1];
			n2tot += ns2 [np - 1];
		} else {
			n2tot = ntop + nbot;
		}
						
		if (ldx == -1) {
			if (id == 0) {
				ldx = ntot;
			} else {
				ldx = ntop + nbot;
			}
		}
		if (ldx < ntop + nbot) {
			throw 0;
		}
		
		matrix_add_scaled (ntop, ntop, 1.0, a, x, lda, ldb);
		matrix_add_scaled (nbot, ntop, 1.0, a + ntop + n, x + ntop, lda, ldb);
		matrix_add_scaled (ntop, nbot, 1.0, a + (ntop + n) * lda, x + ntop * ldb, lda, ldb);
		matrix_add_scaled (nbot, nbot, 1.0, a + (ntop + n) * (lda + 1), x + ntop * (ldb + 1), lda, ldb);
		
		matrix_factorize (n, n, am, ipiv, info, lda);
		
		matrix_solve (n, am, ipiv, al, info, ntop, lda, lda);
		matrix_solve (n, am, ipiv, ar, info, nbot, lda, lda);
		
		matrix_matrix_multiply (ntop, ntop, n, -1.0, at, al, 1.0, x, lda, lda, ldb);
		matrix_matrix_multiply (nbot, ntop, n, -1.0, ab, al, 1.0, x + ntop, lda, lda, ldb);
		matrix_matrix_multiply (ntop, nbot, n, -1.0, at, ar, 1.0, x + ntop * ldx, lda, lda, ldb);
		matrix_matrix_multiply (nbot, nbot, n, -1.0, ab, ar, 1.0, x + ntop * (ldx + 1), lda, lda, ldb);
		
		if (messenger_ptr->get_id () == 0) {
			std::vector <double> buffer (n2tot);
			messenger_ptr->gatherv <double> ((ntop + nbot) * (ntop + nbot), x, &ns2 [0], &buffer [0]);
			for (int i = 0; i < ntot; ++i) {
				scale (ntot, 0.0, x + i * ldx);
			}
			int bcur = 0, xcur = 0;
			for (int i = 0; i < messenger_ptr->get_np () - 1; ++i) {
				for (int j = 0; j < ns [i] + ns [i + 1]; ++j) {
					for (int k = 0; k < ns [i] + ns [i + 1]; ++k) {
						x [xcur + k * ldx + j] += buffer [bcur + k * (ns [i] + ns [i + 1]) + j];
					}
				}
				bcur += (ns [i] + ns [i + 1]) * (ns [i] + ns [i + 1]);
				xcur += ns [i] * (ldx + 1);
			}
			for (int j = 0; j < ns [np - 1]; ++j) {
				for (int k = 0; k < ns [np - 1]; ++k) {
					x [xcur + k * ldx + j] += buffer [bcur + k * ns [np - 1] + j];
				}
			}
			matrix_factorize (ntot, ntot, x, xipiv, info, ldx);
		} else {
			messenger_ptr->gatherv <double> ((ntop + nbot) * (ntop + nbot), x, NULL, NULL);
		}
	}
	
	void p_block_matrix_solve (bases::messenger* messenger_ptr, int n, int ntop, int nbot, double* a, int* ipiv, double* b, double *x, int *xipiv, int *ns, int *info, int nrhs, int lda, int ldx, int ldb) {
		int id = messenger_ptr->get_id (), np = messenger_ptr->get_np ();
		std::vector <int> nsp (np);
		double* am = a + ntop * (lda + 1);
		double* at = a + ntop * lda;
		double* ab = a + ntop * lda + ntop + n;
		double* ar = a + (ntop + n) * lda + ntop;
		double* al = a + ntop;
		
		int ldy = 0;
		if (id == 0) {
			for (int i = 0; i < np - 1; ++i) {
				ldy += ns [i];
				nsp [i] = (ns [i] + ns [i + 1]) * nrhs;
			}
			ldy += ns [np - 1];
			nsp [np - 1] = (ns [np - 1]) * nrhs;
		} else {
			ldy = ntop + nbot;
		}
		
		std::vector <double> y (ldy * nrhs);
		
		if (lda == -1) {
			lda = n + ntop + nbot;
		}
		if (ldb == -1) {
			ldb = n + ntop + nbot;
		}
		if (ldx == -1) {
			ldx = ldy;
		}

		matrix_add_scaled (ntop, nrhs, 1.0, b, &y [0], ldb, ntop + nbot);
		matrix_add_scaled (nbot, nrhs, 1.0, b + ntop + n, &y [ntop], ldb, ntop + nbot);
				
		matrix_solve (n, am, ipiv, b + ntop, info, nrhs, lda, ldb);
		
		matrix_matrix_multiply (ntop, nrhs, n, -1.0, at, b + ntop, 1.0, &y [0], lda, ldb, ntop + nbot);
		matrix_matrix_multiply (nbot, nrhs, n, -1.0, ab, b + ntop, 1.0, &y [ntop], lda, ldb, ntop + nbot);
	
		if (id == 0) {
			int ycur = 0, bcur = 0;
			std::vector <double> buffer (nrhs * 2 * ldy);
			messenger_ptr->gatherv <double> ((ntop + nbot) * nrhs, &y [0], &nsp [0], &buffer [0]);
			utils::scale (nrhs * ldy, 0.0, &y [0]);
			for (int i = 0; i < np; ++i) {
				for (int j = 0; j < nrhs; ++j) {
					for (int k = 0; k < nsp [i] / nrhs; ++k) {
						y [j * ldy + ycur + k] += buffer [bcur + j * nsp [i] / nrhs + k];
					}
				}
				ycur += ns [i];
				bcur += nsp [i];
			}
			
			matrix_solve (ldy, x, xipiv, &y [0], info, nrhs, ldx, ldy);
			
			ycur = 0;
			bcur = 0;
			for (int i = 0; i < np; ++i) {
				for (int j = 0; j < nrhs; ++j) {
					for (int k = 0; k < nsp [i] / nrhs; ++k) {
						buffer [bcur + j * nsp [i] / nrhs + k] = y [j * ldy + ycur + k];
					}
				}
				ycur += ns [i];
				bcur += nsp [i];
			}
			
			messenger_ptr->scatterv <double> ((ntop + nbot) * nrhs, &buffer [0], &nsp [0], &y [0]);
		} else {
			messenger_ptr->gatherv <double> ((ntop + nbot) * nrhs, &y [0], NULL);
			messenger_ptr->scatterv <double> ((ntop + nbot) * nrhs, NULL, &nsp [0], &y [0]);
		}
		
		matrix_copy (ntop, nrhs, &y [0], b, ntop + nbot, ldb);
		matrix_copy (nbot, nrhs, &y [ntop], b + ntop + n, ntop + nbot, ldb);
		
		matrix_matrix_multiply (n, nrhs, ntop, -1.0, al, b, 1.0, b + ntop, lda, ldb, ldb);
		matrix_matrix_multiply (n, nrhs, nbot, -1.0, ar, b + ntop + n, 1.0, b + ntop, lda, ldb, ldb);
	}
	
	void p_block_matrix_factorize (bases::messenger* messenger_ptr, int n, int ntop, int nbot, float* a, int* ipiv, float *x, int *xipiv, int* ns, int *info, int lda, int ldx) {
		int ntot = 0, n2tot = 0, ldb = ntop + nbot, id = messenger_ptr->get_id (), np = messenger_ptr->get_np ();
		float* am = a + ntop * (lda + 1);
		float* at = a + ntop * lda;
		float* ab = a + ntop * lda + ntop + n;
		float* ar = a + (ntop + n) * lda + ntop;
		float* al = a + ntop;
		if (lda == -1) {
			lda = n + ntop + nbot;
		}
				
		std::vector <int> ns2;
		if (id == 0) {
			ns2.resize (np);
			for (int i = 0; i < messenger_ptr->get_np () - 1; ++i) {
				ntot += ns [i];
				ns2 [i] = (ns [i] + ns [i + 1]) * (ns [i] + ns [i + 1]);
				n2tot += ns2 [i];
			}
			ntot += ns [np - 1];
			ns2 [np - 1] = ns [np - 1] * ns [np - 1];
			n2tot += ns2 [np - 1];
		} else {
			n2tot = ntop + nbot;
		}
						
		if (ldx == -1) {
			if (id == 0) {
				ldx = ntot;
			} else {
				ldx = ntop + nbot;
			}
		}
		if (ldx < ntop + nbot) {
			throw 0;
		}
		
		matrix_add_scaled (ntop, ntop, 1.0, a, x, lda, ldb);
		matrix_add_scaled (nbot, ntop, 1.0, a + ntop + n, x + ntop, lda, ldb);
		matrix_add_scaled (ntop, nbot, 1.0, a + (ntop + n) * lda, x + ntop * ldb, lda, ldb);
		matrix_add_scaled (nbot, nbot, 1.0, a + (ntop + n) * (lda + 1), x + ntop * (ldb + 1), lda, ldb);
		
		matrix_factorize (n, n, am, ipiv, info, lda);
		
		matrix_solve (n, am, ipiv, al, info, ntop, lda, lda);
		matrix_solve (n, am, ipiv, ar, info, nbot, lda, lda);
		
		matrix_matrix_multiply (ntop, ntop, n, -1.0, at, al, 1.0, x, lda, lda, ldb);
		matrix_matrix_multiply (nbot, ntop, n, -1.0, ab, al, 1.0, x + ntop, lda, lda, ldb);
		matrix_matrix_multiply (ntop, nbot, n, -1.0, at, ar, 1.0, x + ntop * ldx, lda, lda, ldb);
		matrix_matrix_multiply (nbot, nbot, n, -1.0, ab, ar, 1.0, x + ntop * (ldx + 1), lda, lda, ldb);
		
		if (messenger_ptr->get_id () == 0) {
			std::vector <float> buffer (n2tot);
			messenger_ptr->gatherv <float> ((ntop + nbot) * (ntop + nbot), x, &ns2 [0], &buffer [0]);
			for (int i = 0; i < ntot; ++i) {
				scale (ntot, 0.0, x + i * ldx);
			}
			int bcur = 0, xcur = 0;
			for (int i = 0; i < messenger_ptr->get_np () - 1; ++i) {
				for (int j = 0; j < ns [i] + ns [i + 1]; ++j) {
					for (int k = 0; k < ns [i] + ns [i + 1]; ++k) {
						x [xcur + k * ldx + j] += buffer [bcur + k * (ns [i] + ns [i + 1]) + j];
					}
				}
				bcur += (ns [i] + ns [i + 1]) * (ns [i] + ns [i + 1]);
				xcur += ns [i] * (ldx + 1);
			}
			for (int j = 0; j < ns [np - 1]; ++j) {
				for (int k = 0; k < ns [np - 1]; ++k) {
					x [xcur + k * ldx + j] += buffer [bcur + k * ns [np - 1] + j];
				}
			}
			matrix_factorize (ntot, ntot, x, xipiv, info, ldx);
		} else {
			messenger_ptr->gatherv <float> ((ntop + nbot) * (ntop + nbot), x, NULL, NULL);
		}
	}
	
	void p_block_matrix_solve (bases::messenger* messenger_ptr, int n, int ntop, int nbot, float* a, int* ipiv, float* b, float *x, int *xipiv, int *ns, int *info, int nrhs, int lda, int ldx, int ldb) {
		int id = messenger_ptr->get_id (), np = messenger_ptr->get_np ();
		std::vector <int> nsp (np);
		float* am = a + ntop * (lda + 1);
		float* at = a + ntop * lda;
		float* ab = a + ntop * lda + ntop + n;
		float* ar = a + (ntop + n) * lda + ntop;
		float* al = a + ntop;
		
		int ldy = 0;
		if (id == 0) {
			for (int i = 0; i < np - 1; ++i) {
				ldy += ns [i];
				nsp [i] = (ns [i] + ns [i + 1]) * nrhs;
			}
			ldy += ns [np - 1];
			nsp [np - 1] = (ns [np - 1]) * nrhs;
		} else {
			ldy = ntop + nbot;
		}
		
		std::vector <float> y (ldy * nrhs);
		
		if (lda == -1) {
			lda = n + ntop + nbot;
		}
		if (ldb == -1) {
			ldb = n + ntop + nbot;
		}
		if (ldx == -1) {
			ldx = ldy;
		}

		matrix_add_scaled (ntop, nrhs, 1.0, b, &y [0], ldb, ntop + nbot);
		matrix_add_scaled (nbot, nrhs, 1.0, b + ntop + n, &y [ntop], ldb, ntop + nbot);
				
		matrix_solve (n, am, ipiv, b + ntop, info, nrhs, lda, ldb);
		
		matrix_matrix_multiply (ntop, nrhs, n, -1.0, at, b + ntop, 1.0, &y [0], lda, ldb, ntop + nbot);
		matrix_matrix_multiply (nbot, nrhs, n, -1.0, ab, b + ntop, 1.0, &y [ntop], lda, ldb, ntop + nbot);
	
		if (id == 0) {
			int ycur = 0, bcur = 0;
			std::vector <float> buffer (nrhs * 2 * ldy);
			messenger_ptr->gatherv <float> ((ntop + nbot) * nrhs, &y [0], &nsp [0], &buffer [0]);
			utils::scale (nrhs * ldy, 0.0, &y [0]);
			for (int i = 0; i < np; ++i) {
				for (int j = 0; j < nrhs; ++j) {
					for (int k = 0; k < nsp [i] / nrhs; ++k) {
						y [j * ldy + ycur + k] += buffer [bcur + j * nsp [i] / nrhs + k];
					}
				}
				ycur += ns [i];
				bcur += nsp [i];
			}
			
			matrix_solve (ldy, x, xipiv, &y [0], info, nrhs, ldx, ldy);
			
			ycur = 0;
			bcur = 0;
			for (int i = 0; i < np; ++i) {
				for (int j = 0; j < nrhs; ++j) {
					for (int k = 0; k < nsp [i] / nrhs; ++k) {
						buffer [bcur + j * nsp [i] / nrhs + k] = y [j * ldy + ycur + k];
					}
				}
				ycur += ns [i];
				bcur += nsp [i];
			}
			
			messenger_ptr->scatterv <float> ((ntop + nbot) * nrhs, &buffer [0], &nsp [0], &y [0]);
		} else {
			messenger_ptr->gatherv <float> ((ntop + nbot) * nrhs, &y [0], NULL);
			messenger_ptr->scatterv <float> ((ntop + nbot) * nrhs, NULL, &nsp [0], &y [0]);
		}
		
		matrix_copy (ntop, nrhs, &y [0], b, ntop + nbot, ldb);
		matrix_copy (nbot, nrhs, &y [ntop], b + ntop + n, ntop + nbot, ldb);
		
		matrix_matrix_multiply (n, nrhs, ntop, -1.0, al, b, 1.0, b + ntop, lda, ldb, ldb);
		matrix_matrix_multiply (n, nrhs, nbot, -1.0, ar, b + ntop + n, 1.0, b + ntop, lda, ldb, ldb);
	}
} /* utils */