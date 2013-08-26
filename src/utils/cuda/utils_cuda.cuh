/*!**********************************************************************
 * \file utils_cuda.cuh
 * /Users/justinbrown/Dropbox/spectral_element/src
 * 
 * Created by Justin Brown on 2013-08-26.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef UTILS_CUDA_CUH_AGOKALSC
#define UTILS_CUDA_CUH_AGOKALSC

#define HANDLE_ERROR(status) \
{cudaError_t result = status; \
switch (result) { \
	case cudaErrorMemoryAllocation: FATAL ("Memory Allocation Error."); throw 0; \
	case cudaErrorInvalidValue: FATAL ("Invalid value passed."); throw 0; \
	default: if (status != cudaSuccess) {FATAL ("Other problem."); throw 0;}}}

#define HANDLE_STATUS(status) \
{cublasStatus result = status; \
switch (result) { \
	case CUBLAS_STATUS_NOT_INITIALIZED: printf ("CUBLAS didn't initialize correctly.\n"); throw 0; \
	case CUBLAS_STATUS_ALLOC_FAILED: printf ("CUBLAS allocation failed.\n"); throw 0; \
	case CUBLAS_STATUS_INVALID_VALUE: printf ("CUBLAS unsupported value or parameter.\n"); throw 0; \
	case CUBLAS_STATUS_ARCH_MISMATCH: printf ("CUBLAS feature absent in current architecture.\n"); throw 0; \
	case CUBLAS_STATUS_MAPPING_ERROR: printf ("CUBLAS access to GPU memory failed.\n"); throw 0; \
	case CUBLAS_STATUS_EXECUTION_FAILED: printf ("CUBLAS failed to execute.\n"); throw 0; \
	case CUBLAS_STATUS_INTERNAL_ERROR: printf ("CUBLAS internal operation failed.\n"); throw 0;}}


#endif /* end of include guard: UTILS_CUDA_CUH_AGOKALSC */
