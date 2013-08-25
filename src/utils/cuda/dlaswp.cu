#include <stdio.h>

// Log Change: Adding CUDA functions
const int threadsPerBlock = 256;

__global__ void swap (int n, double* a, int lda, int k1, int k2, int* ipiv, int incx) {
	double temp;
	int ix;
	int ip;
	if (blockIdx.x == 0) {
		if (incx > 0) {
			ix = k1 - 1;
			for (int i = k1 - 1; i < k2; ++i) {
				ip = ipiv [ix] - 1;
				if (ip != i) {
					int j = threadIdx.x;
					if (threadIdx.x == 0) {
					}
					while (j < n) {
						temp = a [i + j * lda];
						a [i + j * lda] = a [ip + j * lda];
						a [ip + j * lda] = temp;
						j += blockDim.x;
					}
					__syncthreads ();
				}
				ix += incx;
			}
		} else {
			ix = (1 - k2) * incx;
			for (int i = k2 - 1; i >= k1 - 1; --i) {
				ip = ipiv [ix] - 1;
				if (ip != i) {
					int j = threadIdx.x;
					if (threadIdx.x == 0) {
					}
					while (j < n) {
						temp = a [i + j * lda];
						a [i + j * lda] = a [ip + j * lda];
						a [ip + j * lda] = temp;
						j += blockDim.x;
					}
					__syncthreads ();
				}
				ix += incx;
			}
		}
	}
}
// End Change

/* Subroutine */ int dlaswpcuda_(int *n, double *a, int *lda, int 
	*k1, int *k2, int *ipiv, int *incx)
{
	// Log Change: Adding CUDA variables
	int threads = min (*n, threadsPerBlock);
	// Original
    /*
    int a_dim1, a_offset, i__1, i__2, i__3, i__4;

    int i__, j, k, i1, i2, n32, ip, ix, ix0, inc;
    double temp;
	*/
	// End Change


/*  -- LAPACK auxiliary routine (version 3.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLASWP performs a series of row interchanges on the matrix A. */
/*  One row interchange is initiated for each of rows K1 through K2 of A. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the matrix of column dimension N to which the row */
/*          interchanges will be applied. */
/*          On exit, the permuted matrix. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A. */

/*  K1      (input) INTEGER */
/*          The first element of IPIV for which a row interchange will */
/*          be done. */

/*  K2      (input) INTEGER */
/*          The last element of IPIV for which a row interchange will */
/*          be done. */

/*  IPIV    (input) INTEGER array, dimension (K2*abs(INCX)) */
/*          The vector of pivot indices.  Only the elements in positions */
/*          K1 through K2 of IPIV are accessed. */
/*          IPIV(K) = L implies rows K and L are to be interchanged. */

/*  INCX    (input) INTEGER */
/*          The increment between successive values of IPIV.  If IPIV */
/*          is negative, the pivots are applied in reverse order. */

/*  Further Details */
/*  =============== */

/*  Modified by */
/*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Interchange row I with row IPIV(I) for each of rows K1 through K2. */

	// Log Change: Adding in CUDA function
	if (incx == 0) {
		return 0;
	}
	
	double* a_dev;
	int* ipiv_dev;
		
	cudaMalloc (&a_dev, sizeof (double) * *n * *lda);
	cudaMalloc (&ipiv_dev, sizeof (int) * *incx * (abs (*k2 - *k1) + 1));
	
	cudaMemcpy (a_dev, a, sizeof (double) * *n * *lda, cudaMemcpyHostToDevice);
	cudaMemcpy (ipiv_dev, ipiv, sizeof (int) * *incx * (abs (*k2 - *k1) + 1), cudaMemcpyHostToDevice);
	
	swap <<<1, threads>>> (*n, a_dev, *lda, *k1, *k2, ipiv_dev, *incx);
	
	cudaDeviceSynchronize ();
	
	cudaMemcpy (a, a_dev, sizeof (double) * *n * *lda, cudaMemcpyDeviceToHost);
	
	cudaFree (a_dev);
	cudaFree (ipiv_dev);
	
	// Original
	/*
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;

    if (*incx > 0) {
		ix0 = *k1;
		i1 = *k1;
		i2 = *k2;
		inc = 1;
	    } else if (*incx < 0) {
		ix0 = (1 - *k2) * *incx + 1;
		i1 = *k2;
		i2 = *k1;
		inc = -1;
    } else {
		return 0;
    }
    * 
    n32 = *n / 32 << 5;
    if (n32 != 0) {
		i__1 = n32;
		for (j = 1; j <= i__1; j += 32) {
		    ix = ix0;
		    i__2 = i2;
		    i__3 = inc;
		    for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
				ip = ipiv[ix];
				if (ip != i__) {
				    i__4 = j + 31;
				    for (k = j; k <= i__4; ++k) {
						temp = a[i__ + k * a_dim1];
						a[i__ + k * a_dim1] = a[ip + k * a_dim1];
						a[ip + k * a_dim1] = temp;
					}
				}
				ix += *incx;
		    }
		}
    }
    if (n32 != *n) {
		++n32;
		ix = ix0;
		i__1 = i2;
		i__3 = inc;
		for (i__ = i1; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
		    ip = ipiv[ix];
		    if (ip != i__) {
				i__2 = *n;
				for (k = n32; k <= i__2; ++k) {
				    temp = a[i__ + k * a_dim1];
				    a[i__ + k * a_dim1] = a[ip + k * a_dim1];
				    a[ip + k * a_dim1] = temp;
				}
		    }
		    ix += *incx;
		}
    }
    */
	// End Changes

    return 0;

/*     End of DLASWP */

} /* dlaswp_ */
