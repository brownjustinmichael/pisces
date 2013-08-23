#include <stdbool.h>
#include <stdio.h>
#include <cuda_runtime.h>

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

// Log change: Changed fortran declarations to C, removed i__1, i__2, i__3

// Log Change: Adding CUDA functions to replace code body
#define HANDLE_ERROR(status) \
{cudaError_t result = status; \
switch (result) { \
	case cudaErrorMemoryAllocation: printf ("Memory Allocation Error.\n"); throw 0; \
	case cudaErrorInvalidValue: printf ("Invalid value passed.\n"); throw 0; \
	default: if (status != cudaSuccess) {printf ("Other problem.\n"); throw 0;}}}

#define TINY 1.e-14

const int threadsPerBlock = 256;

__global__ void zero_matrix (int m, int n, double* a, int lda) {
	int j = blockIdx.x;
	int i = threadIdx.x;
	while (j < n) {
		while (i < m) {
			a [i + j * lda] = 0.0;
			i += blockDim.x;
		}
		j += gridDim.x;
	}
}

__global__ void scale_matrix (int m, int n, double alpha, double* a, int lda) {
	int j = blockIdx.x;
	int i = threadIdx.x;
	while (j < n) {
		while (i < m) {
			a [i + j * lda] *= alpha;
			i += blockDim.x;
		}
		j += gridDim.x;
	}
}

__global__ void operation_left_upper (int m, int n, double* a, int lda, double* b, int ldb, bool nounit) {
	int j = blockIdx.x;
	while (j < n) {
		for (int k = m - 1; k >= 0; --k) {
			int i = threadIdx.x;
			if (i == 0 && nounit && abs (b [k + j * ldb]) != 0.) {
				b [k + j * ldb] /= a [k + k * lda];
			}
			__syncthreads ();
			while (i < k) {
				if (abs (b [k + j * ldb]) != 0.) {
					b [i + j * ldb] -= b [k + j * ldb] * a [i + k * ldb];
				}
				i += blockDim.x;
			}
			__syncthreads ();
		}
		j += gridDim.x;
	}
}

__global__ void operation_left_lower (int m, int n, double* a, int lda, double* b, int ldb, bool nounit) {
	int j = blockIdx.x;
	while (j < n) {
		for (int k = 0; k < m; ++k) {
			int i = threadIdx.x;
			if (i == 0 && nounit && abs (b [k + j * ldb]) != 0.) {
				b [k + j * ldb] /= a [k + k * lda];
			}
			i += k + 1;
			__syncthreads ();
			while (i < m) {
				if (abs (b [k + j * ldb]) != 0.) {
					b [i + j * ldb] -= b [k + j * ldb] * a [i + k * lda];
				}
				i += blockDim.x;
			}
			__syncthreads ();
		}
		j += gridDim.x;
	}
}

__global__ void operation_left_upper_t (int m, int n, double alpha, double* a, int lda, double* b, int ldb, bool nounit) {
	__shared__ double cache [threadsPerBlock];
	int cache_length, half_length;
	int j = blockIdx.x;
	while (j < n) {
		for (int i = 0; i < m; ++i) {
			int k = threadIdx.x;
			if (threadIdx.x == 0) {
				cache [threadIdx.x] = alpha * b [i + j * ldb];
			}
			while (k < i) {
				cache [threadIdx.x] = a [k + i * lda] * b [k + j * ldb];
			}
			k += blockDim.x;
			while (k < i) {
				cache [threadIdx.x] -= a [k + i * lda] * b [k + j * ldb];
				k += blockDim.x;
			}
			__syncthreads ();
			
			k = threadIdx.x;
			cache_length = min (blockDim.x, i + 1);
			while (cache_length > 1) {
				half_length = (cache_length + 1) / 2;
				if (k + half_length < cache_length) {
					cache [k] += cache [k + half_length];
				}
				cache_length = half_length;
				__syncthreads ();
			}
			if (threadIdx.x == 0) {
				if (nounit) {
					b [i + j * ldb] = cache [0] / a [i + i * lda];
				} else {
					b [i + j * ldb] = cache [0];
				}
			}
			__syncthreads ();
		}
		j += gridDim.x;
	}
}

__global__ void operation_left_lower_t (int m, int n, double alpha, double* a, int lda, double* b, int ldb, bool nounit) {
	__shared__ double cache [threadsPerBlock];
	int cache_length, half_length;
	int j = blockIdx.x;
	while (j < n) {
		for (int i = m - 1; i >= 0; --i) {
			int k = threadIdx.x;
			if (threadIdx.x == blockDim.x - 1) {
				cache [threadIdx.x] = alpha * b [i + j * ldb];
			}
			while (k > i) {
				cache [threadIdx.x] = a [k + i * lda] * b [k + j * ldb];
			}
			k += blockDim.x;
			while (k > i) {
				cache [threadIdx.x] -= a [k + i * lda] * b [k + j * ldb];
				k += blockDim.x;
			}
			__syncthreads ();
			
			k = threadIdx.x;
			cache_length = min (blockIdx.x, m - i);
			while (cache_length > 1) {
				half_length = (cache_length + 1) / 2;
				if (k - half_length > blockIdx.x - cache_length) {
					cache [k] += cache [k - half_length];
				}
				cache_length = half_length;
				__syncthreads ();
			}
			if (threadIdx.x == 0) {
				if (nounit) {
					b [i + j * ldb] = cache [blockDim.x - 1] / a [i + i * lda];
				} else {
					b [i + j * ldb] = cache [blockDim.x - 1];
				}
			}
			__syncthreads ();
		}
		j += gridDim.x;
	}
}
// End Change

/* Subroutine */ int dtrsmcuda_(char *side, char *uplo, char *transa, char *diag, 
	int *m, int *n, double *alpha, double *a, int *
	lda, double *b, int *ldb)
{

	/* Local variables */
	static int info;
	static double temp;
	static int i, j, k;
	static bool lside;
	extern bool lsamecuda_(char *, char *);
	static int nrowa;
	static bool upper;
	extern /* Subroutine */ int xerblacuda_(char *, int *);
	static bool nounit;
	
	// Log Change: CUDA variables
	int blocks = min (256, *n), threads = min (threadsPerBlock, *m);
	double* a_dev, *b_dev;
	// End Change


/*  Purpose   
	=======   

	DTRSM  solves one of the matrix equations   

	   op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,   

	where alpha is a scalar, X and B are m by n matrices, A is a unit, or 
  
	non-unit,  upper or lower triangular matrix  and  op( A )  is one  of 
  

	   op( A ) = A   or   op( A ) = A'.   

	The matrix X is overwritten on B.   

	Parameters   
	==========   

	SIDE   - CHARACTER*1.   
			 On entry, SIDE specifies whether op( A ) appears on the left 
  
			 or right of X as follows:   

				SIDE = 'L' or 'l'   op( A )*X = alpha*B.   

				SIDE = 'R' or 'r'   X*op( A ) = alpha*B.   

			 Unchanged on exit.   

	UPLO   - CHARACTER*1.   
			 On entry, UPLO specifies whether the matrix A is an upper or 
  
			 lower triangular matrix as follows:   

				UPLO = 'U' or 'u'   A is an upper triangular matrix.   

				UPLO = 'L' or 'l'   A is a lower triangular matrix.   

			 Unchanged on exit.   

	TRANSA - CHARACTER*1.   
			 On entry, TRANSA specifies the form of op( A ) to be used in 
  
			 the matrix multiplication as follows:   

				TRANSA = 'N' or 'n'   op( A ) = A.   

				TRANSA = 'T' or 't'   op( A ) = A'.   

				TRANSA = 'C' or 'c'   op( A ) = A'.   

			 Unchanged on exit.   

	DIAG   - CHARACTER*1.   
			 On entry, DIAG specifies whether or not A is unit triangular 
  
			 as follows:   

				DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

				DIAG = 'N' or 'n'   A is not assumed to be unit   
									triangular.   

			 Unchanged on exit.   

	M	  - INTEGER.   
			 On entry, M specifies the number of rows of B. M must be at 
  
			 least zero.   
			 Unchanged on exit.   

	N	  - INTEGER.   
			 On entry, N specifies the number of columns of B.  N must be 
  
			 at least zero.   
			 Unchanged on exit.   

	ALPHA  - DOUBLE PRECISION.   
			 On entry,  ALPHA specifies the scalar  alpha. When  alpha is 
  
			 zero then  A is not referenced and  B need not be set before 
  
			 entry.   
			 Unchanged on exit.   

	A	  - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m 
  
			 when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'. 
  
			 Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k 
  
			 upper triangular part of the array  A must contain the upper 
  
			 triangular matrix  and the strictly lower triangular part of 
  
			 A is not referenced.   
			 Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k 
  
			 lower triangular part of the array  A must contain the lower 
  
			 triangular matrix  and the strictly upper triangular part of 
  
			 A is not referenced.   
			 Note that when  DIAG = 'U' or 'u',  the diagonal elements of 
  
			 A  are not referenced either,  but are assumed to be  unity. 
  
			 Unchanged on exit.   

	LDA	- INTEGER.   
			 On entry, LDA specifies the first dimension of A as declared 
  
			 in the calling (sub) program.  When  SIDE = 'L' or 'l'  then 
  
			 LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r' 
  
			 then LDA must be at least max( 1, n ).   
			 Unchanged on exit.   

	B	  - DOUBLE PRECISION array of DIMENSION ( LDB, n ).   
			 Before entry,  the leading  m by n part of the array  B must 
  
			 contain  the  right-hand  side  matrix  B,  and  on exit  is 
  
			 overwritten by the solution matrix  X.   

	LDB	- INTEGER.   
			 On entry, LDB specifies the first dimension of B as declared 
  
			 in  the  calling  (sub)  program.   LDB  must  be  at  least 
  
			 max( 1, m ).   
			 Unchanged on exit.   


	Level 3 Blas routine.   


	-- Written on 8-February-1989.   
	   Jack Dongarra, Argonne National Laboratory.   
	   Iain Duff, AERE Harwell.   
	   Jeremy Du Croz, Numerical Algorithms Group Ltd.   
	   Sven Hammarling, Numerical Algorithms Group Ltd.   



	   Test the input parameters.   

	
   Parameter adjustments   
	   Function Body */

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

	lside = lsamecuda_(side, "L");
	if (lside) {
		nrowa = *m;
	} else {
		nrowa = *n;
	}
	nounit = lsamecuda_(diag, "N");
	upper = lsamecuda_(uplo, "U");

	info = 0;
	if (! lside && ! lsamecuda_(side, "R")) {
		info = 1;
	} else if (! upper && ! lsamecuda_(uplo, "L")) {
		info = 2;
	} else if (! lsamecuda_(transa, "N") && ! lsamecuda_(transa, "T") && ! lsamecuda_(transa, "C")) {
		info = 3;
	} else if (! lsamecuda_(diag, "U") && ! lsamecuda_(diag, "N")) {
		info = 4;
	} else if (*m < 0) {
		info = 5;
	} else if (*n < 0) {
		info = 6;
	} else if (*lda < max(1,nrowa)) {
		info = 9;
	} else if (*ldb < max(1,*m)) {
		info = 11;
	}
	if (info != 0) {
		xerblacuda_("DTRSM ", &info);
		return 0;
	}

/*	 Quick return if possible. */

	if (*n == 0) {
		return 0;
	}

/*	 And when  alpha.eq.zero. */	
		
	if (*alpha == 0.) {
		// Log Change: Adding CUDA function
		HANDLE_ERROR (cudaMalloc (&b_dev, *m * *n * sizeof (double)));
		HANDLE_ERROR (cudaMemcpy (b_dev, b, *m * *n * sizeof (double), cudaMemcpyHostToDevice));

		zero_matrix <<<blocks, threads>>> (*m, *n, b_dev, *ldb);
		
		cudaMemcpy (b, b_dev, *m * *n * sizeof (double), cudaMemcpyDeviceToHost);

		cudaFree (b_dev);
		// Original
		/*
		for (j = 1; j <= *n; ++j) {
			for (i = 1; i <= *m; ++i) {
				B(i,j) = 0.;
			}
		}
		*/
		// End change
	return 0;
	}

/*	 Start the operations. */

	if (lside) {
		if (lsamecuda_(transa, "N")) {

/*		   Form  B := alpha*inv( A )*B. */

			if (upper) {
				// Log Change: Adding CUDA function
				
				HANDLE_ERROR (cudaMalloc (&a_dev, *m * *m * sizeof (double)));
				HANDLE_ERROR (cudaMalloc (&b_dev, *m * *n * sizeof (double)));
	
				HANDLE_ERROR (cudaMemcpy (a_dev, a, *m * *m * sizeof (double), cudaMemcpyHostToDevice));
				HANDLE_ERROR (cudaMemcpy (b_dev, b, *m * *n * sizeof (double), cudaMemcpyHostToDevice));
				
				if (*alpha != 1.0) {
					scale_matrix <<<blocks, threads>>> (*m, *n, *alpha, b_dev, *ldb);
				}
				operation_left_upper <<<blocks, threads>>> (*m, *n, a_dev, *lda, b_dev, *ldb, nounit);
				
				cudaMemcpy (b, b_dev, *m * *n * sizeof (double), cudaMemcpyDeviceToHost);
	
				cudaFree (a_dev);
				cudaFree (b_dev);
				
				// Original
				/*
				for (j = 1; j <= *n; ++j) {
					if (*alpha != 1.) {
						for (i = 1; i <= *m; ++i) {
							B(i,j) = *alpha * B(i,j);
						}
					}
					for (k = *m; k >= 1; --k) {
						if (abs (B(k,j)) > TINY) {
							if (nounit) {
								B(k,j) /= A(k,k);
							}
							for (i = 1; i <= k-1; ++i) {
								B(i,j) -= B(k,j) * A(i,k);
							}
						}
					}
				}
				*/
				// End Change
			} else {
				// Log Change: Adding CUDA function
				double* c;
				c = (double*) malloc (*n * *m);
				HANDLE_ERROR (cudaMalloc (&a_dev, *m * *m * sizeof (double)));
				HANDLE_ERROR (cudaMalloc (&b_dev, *m * *n * sizeof (double)));
	
				HANDLE_ERROR (cudaMemcpy (a_dev, a, *m * *m * sizeof (double), cudaMemcpyHostToDevice));
				HANDLE_ERROR (cudaMemcpy (b_dev, b, *m * *n * sizeof (double), cudaMemcpyHostToDevice));
				
				if (*alpha != 1.0) {
					scale_matrix <<<blocks, threads>>> (*m, *n, *alpha, b_dev, *ldb);
				}
				operation_left_lower <<<blocks, threads>>> (*m, *n, a_dev, *lda, b_dev, *ldb, nounit);
				
				cudaDeviceSynchronize ();
				
				HANDLE_ERROR (cudaMemcpy (b, b_dev, *m * *n * sizeof (double), cudaMemcpyDeviceToHost));

				cudaFree (a_dev);
				cudaFree (b_dev);
				
				cudaDeviceSynchronize ();
				
				// Original
				/*
				for (j = 1; j <= *n; ++j) {
					if (*alpha != 1.) {
						for (i = 1; i <= *m; ++i) {
							B(i,j) = *alpha * B(i,j);
						}
					}
					for (k = 1; k <= *m; ++k) {
						if (abs (B(k,j)) != 0.0) {
							if (nounit) {
								B(k,j) /= A(k,k);
							}
							for (i = k + 1; i <= *m; ++i) {
								printf ("*[%d, %d] %e, %e, %e\n", i-1, k-1, B(i,j),B(k,j), A(i,k));
								B(i,j) -= B(k,j) * A(i,k);
								printf ("*[%d, %d] %e\n", i-1, k-1, B(i,j));
							}
						} else {
							printf ("[%d] Found zero\n", k - 1);
						}
					}
				}
				*/
				// End Change
			}
		} else {
/*		   Form  B := alpha*inv( A' )*B. */
			if (upper) {
				// Log Change: Adding CUDA function, not yet working, kill
				throw 0;
				operation_left_upper_t <<<blocks, threads>>> (*m, *n, *alpha, a_dev, *lda, b_dev, *ldb, nounit);
				// Original
				/*
				for (j = 1; j <= *n; ++j) {
					for (i = 1; i <= *m; ++i) {
						temp = *alpha * B(i,j);
						for (k = 1; k <= i-1; ++k) {
							temp -= A(k,i) * B(k,j);
						}
						if (nounit) {
							temp /= A(i,i);
						}
						B(i,j) = temp;
					}
				}
				*/
				// End Change
			} else {
				// Log Change: Adding CUDA function, not yet working, kill
				throw 0;
				operation_left_lower_t <<<blocks, threads>>> (*m, *n, *alpha, a_dev, *lda, b_dev, *ldb, nounit);
				// Original
				/*
				for (j = 1; j <= *n; ++j) {
					for (i = *m; i >= 1; --i) {
						temp = *alpha * B(i,j);
						for (k = i + 1; k <= *m; ++k) {
							temp -= A(k,i) * B(k,j);
						}
						if (nounit) {
							temp /= A(i,i);
						}
						B(i,j) = temp;
					}
				}
				*/
				// End Change
			}
		}
	} else {
		// Log change: Kill if attempting a right side solve
		throw 0;
		// End change
		if (lsamecuda_(transa, "N")) {

	/*		   Form  B := alpha*B*inv( A ). */

			if (upper) {
			for (j = 1; j <= *n; ++j) {
				if (*alpha != 1.) {
				for (i = 1; i <= *m; ++i) {
					B(i,j) = *alpha * B(i,j);
	/* L170: */
				}
				}
				for (k = 1; k <= j-1; ++k) {
				if (A(k,j) != 0.) {
					for (i = 1; i <= *m; ++i) {
					B(i,j) -= A(k,j) * B(i,k);
	/* L180: */
					}
				}
	/* L190: */
				}
				if (nounit) {
				temp = 1. / A(j,j);
				for (i = 1; i <= *m; ++i) {
					B(i,j) = temp * B(i,j);
	/* L200: */
				}
				}
	/* L210: */
			}
			} else {
			for (j = *n; j >= 1; --j) {
				if (*alpha != 1.) {
				for (i = 1; i <= *m; ++i) {
					B(i,j) = *alpha * B(i,j);
	/* L220: */
				}
				}
				for (k = j + 1; k <= *n; ++k) {
				if (A(k,j) != 0.) {
					for (i = 1; i <= *m; ++i) {
					B(i,j) -= A(k,j) * B(i,k);
	/* L230: */
					}
				}
	/* L240: */
				}
				if (nounit) {
				temp = 1. / A(j,j);
				for (i = 1; i <= *m; ++i) {
					B(i,j) = temp * B(i,j);
	/* L250: */
				}
				}
	/* L260: */
			}
			}
		} else {

	/*		   Form  B := alpha*B*inv( A' ). */

			if (upper) {
			for (k = *n; k >= 1; --k) {
				if (nounit) {
				temp = 1. / A(k,k);
				for (i = 1; i <= *m; ++i) {
					B(i,k) = temp * B(i,k);
	/* L270: */
				}
				}
				for (j = 1; j <= k-1; ++j) {
				if (A(j,k) != 0.) {
					temp = A(j,k);
					for (i = 1; i <= *m; ++i) {
					B(i,j) -= temp * B(i,k);
	/* L280: */
					}
				}
	/* L290: */
				}
				if (*alpha != 1.) {
				for (i = 1; i <= *m; ++i) {
					B(i,k) = *alpha * B(i,k);
	/* L300: */
				}
				}
	/* L310: */
			}
			} else {
			for (k = 1; k <= *n; ++k) {
				if (nounit) {
				temp = 1. / A(k,k);
				for (i = 1; i <= *m; ++i) {
					B(i,k) = temp * B(i,k);
	/* L320: */
				}
				}
				for (j = k + 1; j <= *n; ++j) {
				if (A(j,k) != 0.) {
					temp = A(j,k);
					for (i = 1; i <= *m; ++i) {
					B(i,j) -= temp * B(i,k);
	/* L330: */
					}
				}
	/* L340: */
				}
				if (*alpha != 1.) {
				for (i = 1; i <= *m; ++i) {
					B(i,k) = *alpha * B(i,k);
	/* L350: */
				}
				}
	/* L360: */
			}
			}
		}
	}

	return 0;

/*	 End of DTRSM . */

} /* dtrsmcuda_ */

