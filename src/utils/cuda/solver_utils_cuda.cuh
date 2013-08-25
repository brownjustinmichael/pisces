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
					b [i + j * ldb] -= b [k + j * ldb] * a [i + k * lda];
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
