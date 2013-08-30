#define N 1024
#define M 2
#define S 1

#include "../src/utils/cuda/utils_cublas.cu"
#include "../src/config.cpp"

int main (int argc, char const *argv[])
{
	double* a;
	double* b;
	
	std::vector <cuda::utils::vector <double> > a_devs (S);
	std::vector <cuda::utils::vector <double> > b_devs (S);
	
	cudaEvent_t start, stop;
	float elapsed_time;
	
	std::vector<cudaStream_t> streams (S);
	
	HANDLE_ERROR (cudaEventCreate (&start));
	HANDLE_ERROR (cudaEventCreate (&stop));
	HANDLE_ERROR (cudaEventRecord (start, 0));
	
	for (int i = 0; i < S; ++i) {
		cudaStreamCreate (&streams [i]);
		a_devs [i].resize (N);
		b_devs [i].resize (N);
	}
		
	HANDLE_ERROR (cudaHostAlloc (&a, S * M * N * sizeof (double), cudaHostAllocDefault));
	HANDLE_ERROR (cudaHostAlloc (&b, S * M * N * sizeof (double), cudaHostAllocDefault));
	
	for (int i = 0; i < S * M * N; ++i) {
		a [i] = (double) i;
		b [i] = (double) (i * i);
	}
	
	for (int i = 0; i < S * M * N; ++i) {
		INFO ("In: " << a [i] << " " << b [i]);
	}
	
	for (int i = 0; i < M * S * N; i += S * N) {
		for (int j = 0; j < S; ++j) {
			HANDLE_ERROR (cudaMemcpyAsync (a_devs [j].pointer (), a + j * N + i, N * sizeof (double), cudaMemcpyHostToDevice, streams [j]));
			HANDLE_ERROR (cudaMemcpyAsync (b_devs [j].pointer (), b + j * N + i, N * sizeof (double), cudaMemcpyHostToDevice, streams [j]));
		}

		for (int j = 0; j < S; ++j) {
			cublasSetKernelStream (streams [j]);
			for (int k = 0; k < 20000; ++k) {
				cuda::utils::add_scaled (N, 1.0, a_devs [j].pointer (), b_devs [j].pointer ());
			}
		}

		for (int j = 0; j < S; ++j) {
			HANDLE_ERROR (cudaMemcpyAsync (b + i, b_devs [j].pointer (), N * sizeof (double), cudaMemcpyDeviceToHost, streams [j]));
		}
	}
	
	HANDLE_ERROR (cudaEventRecord (stop, 0));
	
	HANDLE_ERROR (cudaEventSynchronize (stop));
	HANDLE_ERROR (cudaEventElapsedTime (&elapsed_time, start, stop));
	
	for (int i = 0; i < S; ++i) {
		HANDLE_ERROR (cudaStreamSynchronize (streams [i]));
		HANDLE_ERROR (cudaStreamDestroy (streams [i]));
	}
	
	for (int i = 0; i < S * M * N; ++i) {
		INFO ("Out: " << b [i]);
	}
	
	INFO ("Elapsed: " << elapsed_time << " ms.");
	
	HANDLE_ERROR (cudaFreeHost (a));
	HANDLE_ERROR (cudaFreeHost (b));
	
	return 0;
}