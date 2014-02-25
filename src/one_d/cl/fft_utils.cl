#pragma OPENCL EXTENSION cl_khr_fp64 : enable
__kernel void complex_to_real (const int n, __global const double *in, __global double *out) {
	unsigned int dim = 0;
	int tid = 0;
	int multiplier = 1;
	while (dim < get_work_dim ()) {
		tid += multiplier * get_global_id (dim);
		multiplier *= get_global_size (dim);
		++dim;
	}
	while (tid < n) {
		out [tid] = in [tid * 2];
		tid += multiplier;
	}
}

__kernel void symmetrize (const int n, __global const double *in, __global double *out) {
	unsigned int dim = 0;
	int tid = 0;
	int multiplier = 1;
	while (dim < get_work_dim ()) {
		tid += multiplier * get_global_id (dim);
		multiplier *= get_global_size (dim);
		++dim;
	}
	while (tid < n) {
		out [tid] = in [tid];
		out [2 * n - 2 - tid] = in [tid];
		tid += multiplier;
	}
}

__kernel void simple_add (__global const double* A, __global const double* B, __global double* C) {
	C [get_global_id (0)] = A [get_global_id (0)] + B [get_global_id (0)];
}