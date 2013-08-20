#include <cufft.h>
#include <algorithm>
#include <vector>

__global__ void real_to_complex (int n, double* in, cufftDoubleComplex* out) {
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < n) {
		out [tid].x = in [tid];
		out [tid].y = 0.0;
		tid += blockDim.x * gridDim.x;
	}
}

__global__ void complex_to_real (int n, cufftDoubleComplex* in, double* out) {
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < n) {
		out [tid] = in [tid].x;
		tid += blockDim.x * gridDim.x;
	}
}

__global__ void symmetrize (int n, double* data) {
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < n - 1 && tid != 0) {
		data [2 * n - 2 - tid] = data [tid];
		tid += blockDim.x * gridDim.x;
	}
}

class fftw_cosine {
public:
	/*!*******************************************************************
	 * \copydoc bases::transform::transform ()
	 *********************************************************************/
	fftw_cosine (int i_n, double* i_data_in, double* i_data_out) : 
	n (i_n),
	data_in (i_data_in),
	data_out (i_data_out) {
		printf ("Instantiating...");
		if (cudaMalloc ((void **) &data_real, 2 * n * sizeof (double)) != cudaSuccess){
			printf ("Failed to allocate.\n");
			throw 1;	
		}
		if (cudaMalloc ((void **) &data_complex, n * sizeof (cufftDoubleComplex)) != cudaSuccess) {
			printf ("Failed to allocate.\n");
			throw 1;	
		}
		plan = new cufftHandle;
		if (cufftPlan1d(plan, 2 * n - 2, CUFFT_D2Z, 1) != CUFFT_SUCCESS){
			fprintf(stderr, "CUFFT error: Plan creation failed");
			throw 1;	
		}
		printf ("Instantiated.");
	}

	virtual ~fftw_cosine () {
		cufftDestroy (*plan);
		delete plan;
		cudaFree (data_real);
		cudaFree (data_complex);
	}

	/*!*******************************************************************
	 * \copydoc bases::transform::execute ()
	 *********************************************************************/
	void execute ()  {
		for (int i = 0; i < n; ++i) {
			printf ("%d\n", i);
			printf ("thusly %f\n", data_in [i]);
		}

		
		if (cudaMemcpy (data_real, data_in, n * sizeof (double), cudaMemcpyHostToDevice) != cudaSuccess) {
			fprintf (stderr, "FAILURE");
			throw 1;
		}
		
		symmetrize <<<1, std::min (n, 512)>>> (n, data_real);
		
		/* Use the CUFFT plan to transform the signal in place. */
		if (cufftExecD2Z(*plan, data_real, data_complex) != CUFFT_SUCCESS){
			fprintf(stderr, "CUFFT error: ExecC2C Forward failed");
			throw 1;	
		}
		
		if (cudaThreadSynchronize() != cudaSuccess){
			fprintf(stderr, "Cuda error: Failed to synchronize\n");
			throw 1;	
		}
		
		complex_to_real <<<1, std::min (n, 512)>>> (n, data_complex, data_real);
		
		cudaMemcpy (data_out, data_real, n * sizeof (double), cudaMemcpyDeviceToHost);
		
		if (cudaThreadSynchronize() != cudaSuccess){
			fprintf(stderr, "Cuda error: Failed to synchronize\n");
			throw 1;	
		}
		
		for (int i = 0; i < n; ++i) {
			printf ("REAL: %f\n", data_out [i]);
		}
	
	}

private:
	int n;
	double* data_in;
	double* data_out;
	cufftDoubleReal* data_real;
	cufftDoubleComplex* data_complex;
	cufftHandle* plan;
};

int main (int argc, char const *argv[])
{
	int n = 20;
	std::vector <double> data_in (n, 1.0);
	std::vector <double> data_out (n, 0.0);
	
	for (int i = 0; i < n; ++i) {
		data_in [i] = (double) i;
	}
	
	fftw_cosine plan (n, &data_in [0], &data_out [0]);
	plan.execute ();
	return 0;
}