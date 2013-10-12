/*!**********************************************************************
 * \file capability.cu
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-08-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <stdio.h>

int main (int argc, char const *argv[])
{
	if (argc != 2) {
		printf ("Usage: %s major.minor\n", argv [0]);
		return 0;
	}
	
	int major = atoi (&(argv [1] [0]));
	int minor;
	if (argv [1] [1] == '.') {
		minor = atoi (&(argv [1] [2]));
	} else {
		minor = atoi (&(argv [1] [1]));
	}
	
	cudaDeviceProp prop;
	
	int count;
	cudaGetDeviceCount (&count);
	if (cudaGetLastError() != cudaSuccess){
		printf("Cuda error: Failed to get device count.\n");
		throw 1;	
	}
	
	for (int i = 0; i < count; ++i) {
		cudaGetDeviceProperties (&prop, i);
		if (cudaGetLastError() != cudaSuccess){
			printf("Cuda error: Failed to get device properties.\n");
			throw 1;	
		}
		if (prop.major > major) {
			printf ("GPU of compute capability %d.%d can handle CUDA architecture %d.%d.\n", prop.major, prop.minor, major, minor);
			return 0;
		} else if (prop.major < major) {
			printf ("GPU of compute capability %d.%d cannot handle CUDA architecture %d.%d.\n", prop.major, prop.minor, major, minor);
			throw 1;
		} else {
			if (prop.minor >= minor) {
				printf ("GPU of compute capability %d.%d can handle CUDA architecture %d.%d.\n", prop.major, prop.minor, major, minor);
				return 0;
			} else {
				printf ("GPU of compute capability %d.%d cannot handle CUDA architecture %d.%d.\n", prop.major, prop.minor, major, minor);
				throw 1;
			}
		}
	}
	
	return 0;
}