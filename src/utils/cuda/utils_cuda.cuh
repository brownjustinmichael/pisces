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
	case cudaErrorMissingConfiguration: FATAL ("Missing configuration."); throw 0; \
	case cudaErrorMemoryAllocation: FATAL ("Memory Allocation Error."); throw 0; \
	case cudaErrorInitializationError: FATAL ("Initialization Failure."); throw 0; \
	case cudaErrorLaunchFailure: FATAL ("Launch failure."); throw 0; \
	case cudaErrorLaunchTimeout: FATAL ("Launch timeout."); throw 0; \
	case cudaErrorLaunchOutOfResources: FATAL ("Launch out of resources."); throw 0; \
	case cudaErrorInvalidDeviceFunction: FATAL ("Invalid device function."); throw 0; \
	case cudaErrorInvalidConfiguration: FATAL ("Invalid configuration"); throw 0; \
	case cudaErrorInvalidDevice: FATAL ("Invalid device."); throw 0; \
	case cudaErrorInvalidValue: FATAL ("Invalid value passed."); throw 0; \
	case cudaErrorInvalidPitchValue: FATAL ("Invalid pitch value."); throw 0; \
	case cudaErrorInvalidSymbol: FATAL ("Invalid symbol."); throw 0; \
	case cudaErrorMapBufferObjectFailed: FATAL ("Map buffer object failed."); throw 0; \
	case cudaErrorUnmapBufferObjectFailed: FATAL ("Unmap buffer object failed."); throw 0; \
	case cudaErrorInvalidHostPointer: FATAL ("Invalid host pointer."); throw 0; \
	case cudaErrorInvalidDevicePointer: FATAL ("Invalid device pointer."); throw 0; \
	case cudaErrorInvalidTexture: FATAL ("Invalid texture."); throw 0; \
	case cudaErrorInvalidTextureBinding: FATAL ("Invalid texture binding."); throw 0; \
	case cudaErrorInvalidChannelDescriptor: FATAL ("Invalid channel descriptor."); throw 0; \
	case cudaErrorInvalidMemcpyDirection: FATAL ("Invalid Memcpy direction."); throw 0; \
	case cudaErrorInvalidFilterSetting: FATAL ("Invalid filter setting."); throw 0; \
	case cudaErrorInvalidNormSetting: FATAL ("Invalid norm setting."); throw 0; \
	case cudaErrorCudartUnloading: FATAL ("Cudart unloading error"); throw 0; \
	case cudaErrorUnknown: FATAL ("Unknown error."); throw 0; \
	default: if (result != cudaSuccess) {FATAL ("Other problem: " << result << " " << cudaErrorInvalidDevicePointer); throw 0;}}}

#define HANDLE_STATUS(status) \
{cublasStatus_t result = status; \
switch (result) { \
	case CUBLAS_STATUS_NOT_INITIALIZED: FATAL ("CUBLAS didn't initialize correctly."); throw 0; \
	case CUBLAS_STATUS_ALLOC_FAILED: FATAL ("CUBLAS allocation failed."); throw 0; \
	case CUBLAS_STATUS_INVALID_VALUE: FATAL ("CUBLAS unsupported value or parameter."); throw 0; \
	case CUBLAS_STATUS_ARCH_MISMATCH: FATAL ("CUBLAS feature absent in current architecture."); throw 0; \
	case CUBLAS_STATUS_MAPPING_ERROR: FATAL ("CUBLAS access to GPU memory failed."); throw 0; \
	case CUBLAS_STATUS_EXECUTION_FAILED: FATAL ("CUBLAS failed to execute."); throw 0; \
	case CUBLAS_STATUS_INTERNAL_ERROR: FATAL ("CUBLAS internal operation failed."); throw 0;}}


#endif /* end of include guard: UTILS_CUDA_CUH_AGOKALSC */
