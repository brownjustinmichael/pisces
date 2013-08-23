/*
 * Copyright 1993-2013 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */


// System includes
#include <stdio.h>

__global__ void testKernel(int val)
{
	printf("value %d\n", val);
    
}

int main(int argc, char **argv)
{
    printf("printf() is called. Output:\n\n");
	testKernel<<<1, 1>>>(10);

    cudaDeviceSynchronize();

    return EXIT_SUCCESS;
}

