/*!**********************************************************************
 * \file linalg.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-08-07.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "linalg.hpp"
#include "exceptions.hpp"
#include <omp.h>

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
 * \param a A float matrix to be overwritten with its LU decomposition
 * \param lda A pointer to the integer number of leading dimension of a
 * \param ipiv An integer array to contain the pivot indices
 * \param b The float right hand side array, overwritten with solution
 * \param ldb A pointer to the integer number of leading dimension of b
 * \param info A pointer to an integer indicating success (0 for successful exit)
 *********************************************************************/
extern "C" void sgetrs_ (char *trans, int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, int *info);

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

/*!**********************************************************************
 * \brief Function from LAPACK that factorizes a banded matrix
 * 
 * \param m A pointer to the number of rows in the matrix
 * \param n A pointer to the number of columns in the matrix
 * \param kl A pointer to the number of supdiagonal rows in the matrix
 * \param ku A pointer to the number of superdiagonal rows in the matrix
 * \param ab A pointer to the double matrix (must have dimensions n * (2 * kl + ku + 1))
 * \param ldab A pointer to the integer leading dimension of ab
 * \param ipiv An integer array to contain the pivot indices
 * \param info A pointer to an integer indicating success (0 for successful exit)
 ************************************************************************/
extern "C" void dgbtrf_ (int *m, int *n, int *kl, int *ku, double *ab, int *ldab, int *ipiv, int *info);

/*!**********************************************************************
 * \brief Function from LAPACK that factorizes a banded matrix
 * 
 * \param trans A pointer to transposition character ("N" for not transposed, "T" for transposed)
 * \param n A pointer to the number of columns in the matrix
 * \param kl A pointer to the number of supdiagonal rows in the matrix
 * \param ku A pointer to the number of superdiagonal rows in the matrix
 * \param nrhs A pointer to the number of right hand sides of the equation
 * \param ab A pointer to the double matrix (must have dimensions n * (2 * kl + ku + 1))
 * \param ldab A pointer to the integer leading dimension of ab
 * \param ipiv An integer array to contain the pivot indices
 * \param b A pointer to the double array right hand side
 * \param ldb A pointer to the integer leading dimension of b
 * \param info A pointer to an integer indicating success (0 for successful exit)
 ************************************************************************/
extern "C" void dgbtrs_ (char *trans, int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info);

/*!**********************************************************************
 * \brief Function from LAPACK that factorizes a tridiagonal matrix
 * 
 * \param n A pointer to the number of columns in the matrix
 * \param dl A pointer to the float subdiagonal elements of the matrix
 * \param d A pointer to the float diagonal elements of the matrix
 * \param du A pointer to the float superdiagonal elements of the matrix
 * \param du2 A pointer to the float supersuperdiagonal elements of the matrix
 * \param ipiv An integer array to contain the pivot indices
 * \param info A pointer to an integer indicating success (0 for successful exit)
 ************************************************************************/
extern "C" void sgttrf_ (int *n, float *dl, float *d, float *du, float *du2, int *ipiv, int *info);

/*!**********************************************************************
 * \brief Function from LAPACK that factorizes a tridiagonal matrix
 * 
 * \param n A pointer to the number of columns in the matrix
 * \param dl A pointer to the double subdiagonal elements of the matrix
 * \param d A pointer to the double diagonal elements of the matrix
 * \param du A pointer to the double superdiagonal elements of the matrix
 * \param du2 A pointer to the double supersuperdiagonal elements of the matrix
 * \param ipiv An integer array to contain the pivot indices
 * \param info A pointer to an integer indicating success (0 for successful exit)
 ************************************************************************/
extern "C" void dgttrf_ (int *n, double *dl, double *d, double *du, double *du2, int *ipiv, int *info);

/*!**********************************************************************
 * \brief Function from LAPACK that solves a tridiagonal matrix
 * 
 * \param trans A pointer to transposition character ("N" for not transposed, "T" for transposed)
 * \param n A pointer to the number of columns in the matrix
 * \param nrhs A pointer to the integer number of right hand sides in b
 * \param dl A pointer to the float subdiagonal elements of the matrix
 * \param d A pointer to the float diagonal elements of the matrix
 * \param du A pointer to the float superdiagonal elements of the matrix
 * \param du2 A pointer to the float supersuperdiagonal elements of the matrix
 * \param ipiv An integer array to contain the pivot indices
 * \param b A pointer to the float array right hand side
 * \param ldb A pointer to the integer leading dimension of b
 * \param info A pointer to an integer indicating success (0 for successful exit)
 ************************************************************************/
extern "C" void sgttrs_ (char *trans, int *n, int *nrhs, float *dl, float *d, float *du, float *du2, int *ipiv, float *b, int *ldb, int *info);

/*!**********************************************************************
 * \brief Function from LAPACK that solves a tridiagonal matrix
 * 
 * \param trans A pointer to transposition character ("N" for not transposed, "T" for transposed)
 * \param n A pointer to the number of columns in the matrix
 * \param nrhs A pointer to the integer number of right hand sides in b
 * \param dl A pointer to the double subdiagonal elements of the matrix
 * \param d A pointer to the double diagonal elements of the matrix
 * \param du A pointer to the double superdiagonal elements of the matrix
 * \param du2 A pointer to the double supersuperdiagonal elements of the matrix
 * \param ipiv An integer array to contain the pivot indices
 * \param b A pointer to the float array right hand side
 * \param ldb A pointer to the integer leading dimension of b
 * \param info A pointer to an integer indicating success (0 for successful exit)
 ************************************************************************/
extern "C" void dgttrs_ (char *trans, int *n, int *nrhs, double *dl, double *d, double *du, double *du2, int *ipiv, double *b, int *ldb, int *info);

/*!**********************************************************************
 * \brief Function from LAPACK that factorizes and solves a tridiagonal matrix
 * 
 * \param n A pointer to the number of columns in the matrix
 * \param nrhs A pointer to the integer number of right hand sides in b
 * \param sub A pointer to the float subdiagonal elements of the matrix (overwritten on output)
 * \param diag A pointer to the float diagonal elements of the matrix (overwritten on output)
 * \param sup A pointer to the float superdiagonal elements of the matrix (overwritten on output)
 * \param b A pointer to the float array right hand side
 * \param ldb A pointer to the integer leading dimension of b
 * \param info A pointer to an integer indicating success (0 for successful exit)
 * 
 * Warning: the matrix will be overwritted during output and will need to be reset before this function is run again
 ************************************************************************/
extern "C" void sgtsv_ (int *n, int *nrhs, float *sub, float *diag, float *sup, float *b, int *ldb, int *info);

/*!**********************************************************************
 * \brief Function from LAPACK that factorizes and solves a tridiagonal matrix
 * 
 * \param n A pointer to the number of columns in the matrix
 * \param nrhs A pointer to the integer number of right hand sides in b
 * \param sub A pointer to the double subdiagonal elements of the matrix (overwritten on output)
 * \param diag A pointer to the double diagonal elements of the matrix (overwritten on output)
 * \param sup A pointer to the double superdiagonal elements of the matrix (overwritten on output)
 * \param b A pointer to the double array right hand side
 * \param ldb A pointer to the integer leading dimension of b
 * \param info A pointer to an integer indicating success (0 for successful exit)
 * 
 * Warning: the matrix will be overwritted during output and will need to be reset before this function is run again
 ************************************************************************/
extern "C" void dgtsv_ (int *n, int *nrhs, double *sub, double *diag, double *sup, double *b, int *ldb, int *info);

namespace linalg
{
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
			FATAL ("ERROR " << *info);
			throw exceptions::cannot_factor ();
		}
	}
	
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
			FATAL ("ERROR " << *info);
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
				
		int threads = omp_get_max_threads ();
		#pragma omp parallel for
		for (int i = 0; i < threads; ++i) {
			int tnrhs = nrhs / threads + (i < nrhs % threads ? 1 : 0);
			int tldb = threads * ldb;
			sgetrs_ (&charN, &n, &tnrhs, a, &lda, ipiv, b + i * ldb, &tldb, info);
		}
		
		if (*info != 0) {
			throw exceptions::cannot_solve ();
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
		
		int threads = omp_get_max_threads ();
		#pragma omp parallel for
		for (int i = 0; i < threads; ++i) {
			int tnrhs = nrhs / threads + (i < nrhs % threads ? 1 : 0);
			int tldb = threads * ldb;
			dgetrs_ (&charN, &n, &tnrhs, a, &lda, ipiv, b + i * ldb, &tldb, info);
		}
		
		if (*info != 0) {
			throw exceptions::cannot_solve ();
		}
	}
	
	void matrix_banded_factorize (int m, int n, int kl, int ku, double* a, int *ipiv, int *info, int lda) {
		int iinfo;
		if (!info) {
			info = &iinfo;
		}
		
		if (m == 0 || n == 0) {
			return;
		}
		
		if (lda == -1) {
			lda = 2 * kl + ku + 1;
		}
		
		dgbtrf_ (&m, &n, &kl, &ku, a, &lda, ipiv, info);
			
		if (*info != 0) {
			FATAL ("ERROR " << *info);
			throw exceptions::cannot_factor ();
		}
	}
	
	void matrix_banded_solve (int n, int kl, int ku, double* a, int* ipiv, double* b, int *info, int nrhs, int lda, int ldb) {
		char charN = 'N';
		int iinfo;
		
		if (!info) {
			info = &iinfo;
		}
		
		if (n == 0) {
			return;
		}
		
		if (lda == -1) {
			lda = 2 * kl + ku + 1;
		}
		
		if (ldb == -1) {
			ldb = n;
		}
		
		dgbtrs_ (&charN, &n, &kl, &ku, &nrhs, a, &lda, ipiv, b, &ldb, info);
		
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
	
	void tridiagonal_factorize (int n, double *dl, double *d, double *du, double *du2, int *ipiv, int *info) {
		int iinfo;
		if (!info) {
			info = &iinfo;
		}
				
		dgttrf_ (&n, dl, d, du, du2, ipiv, info);
		
		if (*info != 0 && *info != n) {
			FATAL ("Info is " << *info);
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
} /* linalg */