#include <stdbool.h>
#include <stdio.h>
#include "solver_utils_cuda.cuh"

int xerblacuda_(char *, int *);

// Log Change: Rewriting names and variables
int solve_lu (int n, int nrhs, double *a, int lda, int* ipiv, double*b, int ldb, int* info) {
// Original
/* 
static int c__1 = 1;
static double c_b12 = 1.;
static int c_n1 = -1;

* int dgetrs_(char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info)
{
int a_dim1, a_offset, b_dim1, b_offset, alt_info;

	extern int dlaswp_(int *n, double *a, int *lda, int *k1, int *k2, int *ipiv, int *incx);
    extern bool lsame_(char *, char *);
    extern int dtrsm_(char *, char *, char *, char *, int *, int *, double *, double *, int *, double *, int *)
    xerbla_(char *, int *)
    dlaswp_(int *, double *, int *, int *, int *, int *, int *);
    bool notran;
	*/
// End Change


/*  -- LAPACK routine (version 3.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGETRS solves a system of linear equations */
/*     A * X = B  or  A' * X = B */
/*  with a general N-by-N matrix A using the LU factorization computed */
/*  by DGETRF. */

/*  Arguments */
/*  ========= */

/*  TRANS   (input) CHARACTER*1 */
/*          Specifies the form of the system of equations: */
/*          = 'N':  A * X = B  (No transpose) */
/*          = 'T':  A'* X = B  (Transpose) */
/*          = 'C':  A'* X = B  (Conjugate transpose = Transpose) */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  NRHS    (input) INTEGER */
/*          The number of right hand sides, i.e., the number of columns */
/*          of the matrix B.  NRHS >= 0. */

/*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*          The factors L and U from the factorization A = P*L*U */
/*          as computed by DGETRF. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,N). */

/*  IPIV    (input) INTEGER array, dimension (N) */
/*          The pivot indices from DGETRF; for 1<=i<=N, row i of the */
/*          matrix was interchanged with row IPIV(i). */

/*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS) */
/*          On entry, the right hand side matrix B. */
/*          On exit, the solution matrix X. */

/*  LDB     (input) INTEGER */
/*          The leading dimension of the array B.  LDB >= max(1,N). */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

/*     Quick return if possible */
    if (n == 0 || nrhs == 0) {
	return 0;
    }
	
	// Log Change: Only needed for No Transpose, replaced cody body with CUDA calls
    *info = 0;
	int alt_info;
    if (n < 0) {
		*info = -2;
    } else if (nrhs < 0) {
		*info = -3;
    } else if (lda < max(1,n)) {
		*info = -5;
    } else if (ldb < max(1,n)) {
		*info = -8;
    }
    if (*info != 0) {
		alt_info = -(*info);
		xerblacuda_("DGETRS", &alt_info);
		return 0;
	}
	
	swap <<<1, min (nrhs, 256)>>> (nrhs, b, ldb, 1, n, ipiv, 1);
	operation_left_lower <<<min (nrhs, 1024), min (n, 256)>>> (n, nrhs, a, lda, b, ldb, false);
	operation_left_upper <<<min (nrhs, 1024), min (n, 256)>>> (n, nrhs, a, lda, b, ldb, true);
	
	
// 	double *a_dev, *b_dev;
// 	int *ipiv_dev;
// 	
// 	cudaMalloc (&a_dev, sizeof (double) * n * lda);
// 	cudaMalloc (&b_dev, sizeof (double) * nrhs * ldb);
// 	cudaMalloc (&ipiv_dev, sizeof (int) * n);
// 	
// 	cudaMemcpy (a_dev, a, sizeof (double) * n * lda, cudaMemcpyHostToDevice);
// 	cudaMemcpy (b_dev, b, sizeof (double) * nrhs * ldb, cudaMemcpyHostToDevice);
// 	cudaMemcpy (ipiv_dev, ipiv, sizeof (int) * n, cudaMemcpyHostToDevice);
// 	
// 	swap <<<1, min (nrhs, 256)>>> (nrhs, b_dev, ldb, 1, n, ipiv_dev, 1);
// 	
// /*        Solve L*X = B, overwriting B with X. */
// 
// 	operation_left_lower <<<min (nrhs, 1024), min (n, 256)>>> (n, nrhs, a_dev, lda, b_dev, ldb, false);
// 	
// /*        Solve U*X = B, overwriting B with X. */
// 
// 	operation_left_upper <<<min (nrhs, 1024), min (n, 256)>>> (n, nrhs, a_dev, lda, b_dev, ldb, true);
// 		
// 	cudaMemcpy (b, b_dev, sizeof (double) * nrhs * ldb, cudaMemcpyDeviceToHost);
// 	
// 	cudaDeviceSynchronize ();
// 	
// 	cudaFree (ipiv_dev);
// 	cudaFree (a_dev);
// 	cudaFree (b_dev);
	// Original
	/*
    // a_dim1 = *lda;
    // a_offset = 1 + a_dim1;
	// a -= a_offset;
	// --ipiv;
    // b_dim1 = *ldb;
    // b_offset = 1 + b_dim1;
	// b -= b_offset;
	
    *info = 0;
    notran = lsamecuda_(trans, "N");
    if (! notran && ! lsamecuda_(trans, "T") && ! lsamecuda_(
	    trans, "C")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	alt_info = -(*info);
	xerblacuda_("DGETRS", &alt_info);
	return 0;
    }

    if (notran) {
		
		dlaswpcuda_(nrhs, b, ldb, &c__1, n, ipiv, &c__1);
		
		dtrsmcuda_("Left", "Lower", "No transpose", "Unit", n, nrhs, &c_b12, a, lda, b, ldb);

		dtrsmcuda_("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b12, a, lda, b, ldb);
	
    } else {

	dtrsmcuda_("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b12, &a[
		a_offset], lda, &b[b_offset], ldb);

	dtrsmcuda_("Left", "Lower", "Transpose", "Unit", n, nrhs, &c_b12, &a[
		a_offset], lda, &b[b_offset], ldb);

	dlaswpcuda_(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c_n1);
    }
    */
	// End Change

    return 0;

/*     End of DGETRS */

} /* dgetrs_ */
