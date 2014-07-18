/*!**********************************************************************
 * \file solver_lapack.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "solver_utils.hpp"
#include "exceptions.hpp"

/*!*******************************************************************
 * \brief Function from LAPACK that factorizes the matrix a by LU decomposition
 * 
 * \param m A pointer to the number of rows in a
 * \param n A pointer to the number of columns in a
 * \param a A double matrix to be overwritten with its LU decomposition
 * \param lda A double pointer to the integer leading dimension of a
 * \param ipiv An integer array to contain the pivot indices
 * \param info A pointer to an integer indicating success (0 for successful exit)
 *********************************************************************/
extern "C" void dgetrf_ (int *m, int *n, double *a, int *lda, int *ipiv, int *info);

/*!*******************************************************************
 * \brief Function from LAPACK that solves a factorized matrix equation
 * 
 * \param trans A pointer to transposition character ("N" for not transposed, "T" for transposed)
 * \param n A pointer to the number of columns in a
 * \param nrhs A pointer to the number of right hand sides
 * \param a A double matrix to be overwritten with its LU decomposition
 * \param lda A pointer to the integer number of leading dimension of a
 * \param ipiv An integer array to contain the pivot indices
 * \param b The double right hand side array, overwritten with solution
 * \param ldb A pointer to the integer number of leading dimension of b
 * \param info A pointer to an integer indicating success (0 for successful exit)
 *********************************************************************/
extern "C" void dgetrs_ (char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

/*!*******************************************************************
 * \brief Function from LAPACK that factorizes the matrix a by LU decomposition
 * 
 * \param m A pointer to the number of rows in a
 * \param n A pointer to the number of columns in a
 * \param a A float matrix to be overwritten with its LU decomposition
 * \param lda A float pointer to the integer leading dimension of a
 * \param ipiv An integer array to contain the pivot indices
 * \param info A pointer to an integer indicating success (0 for successful exit)
 *********************************************************************/
extern "C" void sgetrf_ (int *m, int *n, float *a, int *lda, int *ipiv, int *info);

/*!*******************************************************************
 * \brief Function from LAPACK that solves a factorized matrix equation
 * 
 * \param trans A pointer to transposition character ("N" for not transposed, "T" for transposed)
 * \param n A pointer to the number of columns in a
 * \param nrhs A pointer to the number of right hand sides
 * \param a A float matrix to be overwritten with its LU decomposition
 * \param lda A pointer to the integer number of leading dimension of a
 * \param ipiv An integer array to contain the pivot indices
 * \param b The float right hand side array, overwritten with solution
 * \param ldb A pointer to the integer number of leading dimension of b
 * \param info A pointer to an integer indicating success (0 for successful exit)
 *********************************************************************/
extern "C" void sgetrs_ (char *trans, int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, int *info);

extern "C" void sgetrs_ (char *trans, int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, int *info);

extern "C" void sgtsv_ (int *n, int *nrhs, float *sub, float *diag, float *sup, float *b, int *ldb, int *info);

extern "C" void dgtsv_ (int *n, int *nrhs, double *sub, double *diag, double *sup, double *b, int *ldb, int *info);

extern "C" void dgttrf_ (int *n, double *dl, double *d, double *du, double *du2, int *ipiv, int *info);

extern "C" void dgttrs_ (char *trans, int *n, int *nrhs, double *dl, double *d, double *du, double *du2, int *ipiv, double *b, int *ldb, int *info);

extern "C" void sgttrf_ (int *n, float *dl, float *d, float *du, float *du2, int *ipiv, int *info);

extern "C" void sgttrs_ (char *trans, int *n, int *nrhs, float *dl, float *d, float *du, float *du2, int *ipiv, float *b, int *ldb, int *info);


namespace utils
{
	void matrix_factorize (int m, int n, double* a, int *ipiv, int *info, int lda) {
		int iinfo;
		if (!info) {
			info = &iinfo;
		}
		
		if (m == 0 || n == 0) {
			return;
		}
		
		if (lda == -1) {
			lda = m;
		}
		
		dgetrf_ (&m, &n, a, &lda, ipiv, info);
		
		if (*info != 0) {
			throw exceptions::cannot_factor ();
		}
	}
	
	void matrix_solve (int n, double* a, int* ipiv, double* b, int *info, int nrhs, int lda, int ldb) {
		char charN = 'N';
		int iinfo;
		
		if (!info) {
			info = &iinfo;
		}
		
		if (n == 0) {
			return;
		}
		
		if (lda == -1) {
			lda = n;
		}
		
		if (ldb == -1) {
			ldb = n;
		}
		
		dgetrs_ (&charN, &n, &nrhs, a, &lda, ipiv, b, &ldb, info);
		
		if (*info != 0) {
			throw exceptions::cannot_solve ();
		}
	}
	
	void matrix_factorize (int m, int n, float* a, int *ipiv, int *info, int lda) {
		int iinfo;
		if (!info) {
			info = &iinfo;
		}
		
		if (m == 0 || n == 0) {
			return;
		}
		
		if (lda == -1) {
			lda = m;
		}
		
		sgetrf_ (&m, &n, a, &lda, ipiv, info);
		
		if (*info != 0) {
			throw exceptions::cannot_factor ();
		}
	}
	
	void matrix_solve (int n, float* a, int* ipiv, float* b, int *info, int nrhs, int lda, int ldb) {
		char charN = 'N';
		int iinfo;
		
		if (!info) {
			info = &iinfo;
		}
		
		if (n == 0 || nrhs == 0) {
			return;
		}
		
		if (lda == -1) {
			lda = n;
		}
		
		if (ldb == -1) {
			ldb = n;
		}
				
		sgetrs_ (&charN, &n, &nrhs, a, &lda, ipiv, b, &ldb, info);
		
		if (*info != 0) {
			throw exceptions::cannot_solve ();
		}
	}
	
	void tridiagonal_factorize (int n, double *dl, double *d, double *du, double *du2, int *ipiv, int *info) {
		int iinfo;
		if (!info) {
			info = &iinfo;
		}
		
		// printf ("POINTER: %i %p %p %p %p %p\n", n, dl, d, du, du2, info);
		
		// for (int i = 0; i < n; ++i) {
		// 	printf ("%f %f %f \n", dl [i], d [i], du [i]);
		// }		
		dgttrf_ (&n, dl, d, du, du2, ipiv, info);
		
		if (*info != 0) {
			FATAL ("Info is " << *info);
			throw exceptions::cannot_factor ();
		}
	}
	
	void tridiagonal_solve (int n, double *dl, double *d, double *du, double *du2, int *ipiv, double *b, int *info, int nrhs, int ldb) {
		char charN = 'N';
		int iinfo;
		
		if (n == 0) {
			return;
		}
		
		if (!info) {
			info = &iinfo;
		}
		if (ldb == -1) {
			ldb = n;
		}

		dgttrs_ (&charN, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, info);
		
		if (*info != 0) {
			throw exceptions::cannot_solve ();
		}
	}

	void tridiagonal_factorize (int n, float *dl, float *d, float *du, float *du2, int *ipiv, int *info) {
		int iinfo;
		
		if (!info) {
			info = &iinfo;
		}
		sgttrf_ (&n, dl, d, du, du2, ipiv, info);
		
		if (*info != 0) {
			throw exceptions::cannot_factor ();
		}
	}
	
	void tridiagonal_solve (int n, float *dl, float *d, float *du, float *du2, int *ipiv, float *b, int *info, int nrhs, int ldb) {
		char charN = 'N';
		int iinfo;
		
		if (n == 0) {
			return;
		}
		
		if (!info) {
			info = &iinfo;
		}
		if (ldb == -1) {
			ldb = n;
		}
		
		sgttrs_ (&charN, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, info);
		
		if (*info != 0) {
			throw exceptions::cannot_solve ();
		}
	}

	void tridiagonal_solve (int n, float *sub, float *diag, float *sup, float *b, int *info, int nrhs, int ldb) {
		int iinfo;
		if (!info) {
			info = &iinfo;
		}
		
		if (ldb == -1) {
			ldb = n;
		}
	
		sgtsv_ (&n, &nrhs, sub, diag, sup, b, &ldb, info);
		
		if (*info != 0) {
			throw exceptions::cannot_solve ();
		}
	}
	
	void tridiagonal_solve (int n, double *sub, double *diag, double *sup, double *b, int *info, int nrhs, int ldb) {
		int iinfo;
		if (!info) {
			info = &iinfo;
		}
		
		if (ldb == -1) {
			ldb = n;
		}
	
		dgtsv_ (&n, &nrhs, sub, diag, sup, b, &ldb, info);
		
		if (*info != 0) {
			throw exceptions::cannot_solve ();
		}
	}
	
	void tridiagonal_direct_solve (int n, double *dl, double *d, double *du, double* du2, double *b, int nrhs, int ldb) {
		if (n == 0) {
			return;
		}
		
		if (ldb == -1) {
			ldb = n;
		}

		dl -= 1;
		
		copy (n, du, du2);
		
		if (d [0] == 0.0) {
			FATAL ("Algorithm can't handle 0 in first element of diagonal.");
			throw 0;
		}
		
		du2 [0] /= d [0];
		for (int i = 0; i < nrhs; ++i) {
			b [i * ldb] /= d [0];
		}
		for (int j = 1; j < n; ++j) {
			du2 [j] /= d [j] - dl [j] * du2 [j - 1];
			for (int i = 0; i < nrhs; ++i) {
				b [i * ldb + j] = (b [i * ldb + j] - dl [j] * b [i * ldb + j - 1]) / (d [j] - dl [j] * du2 [j - 1]);
				if (b [i * ldb + j] != b [i * ldb + j]) {
					DEBUG ("FOUND NAN " << b [i * ldb + j] << " " << dl [j] << " " << b [i * ldb + j - 1] << " " << d [j] << " " << du2 [j - 1]);
					throw 0;
				}
			}
		}
		
		for (int j = n - 2; j >= 0; --j) {
			for (int i = 0; i < nrhs; ++i) {
				b [i * ldb + j] -= du2 [j] * b [i * ldb + j + 1];
			}
		}
	}

	void tridiagonal_direct_solve (int n, float *dl, float *d, float *du, float* du2, float *b, int nrhs, int ldb) {
		if (n == 0) {
			return;
		}
		
		if (ldb == -1) {
			ldb = n;
		}
		
		dl -= 1;
		copy (n, du, du2);
		
		if (d [0] == 0.0) {
			FATAL ("Algorithm can't handle 0 in first element of diagonal.");
			throw 0;
		}
		du2 [0] /= d [0];
		for (int i = 0; i < nrhs; ++i) {
			b [i * ldb] /= d [0];
		}
		for (int j = 1; j < n; ++j) {
			du2 [j] /= d [j] - dl [j] * du2 [j - 1];
			for (int i = 0; i < nrhs; ++i) {
				b [i * ldb + j] = (b [i * ldb + j] - dl [j] * b [i * ldb + j - 1]) / (d [j] - dl [j] * du2 [j - 1]);
			}
		}
		
		for (int j = n - 2; j >= 0; ++j) {
			for (int i = 0; i < nrhs; ++i) {
				b [i * ldb + j] -= du2 [j] * b [i * ldb + j + 1];
				if (b [i * ldb + j] != b [i * ldb + j]) {
					FATAL ("FOUND NAN IN TRIDIAG DIRECT SOLVE");
					throw 0;
				}
			}
		}
	}

	void diagonal_solve (int n, float *a, float *b, int inca, int incb) {
		for (int i = 0; i < n; ++i) {
			b [i * incb] /= a [i * inca]; 
		}
	}

	void diagonal_solve (int n, double *a, double *b, int inca, int incb) {
		for (int i = 0; i < n; ++i) {
			b [i * incb] /= a [i * inca]; 
		}
	}
} /* utils */