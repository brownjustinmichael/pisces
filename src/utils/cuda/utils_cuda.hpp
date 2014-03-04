/*!**********************************************************************
 * \file utils_cuda.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-26.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef UTILS_CUDA_CUH_AGOKALSC
#define UTILS_CUDA_CUH_AGOKALSC

#include <stddef.h>

/*
	TODO HANDLE_ERROR should be in a cuh file
*/

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
	
#define CUBLAS_HANDLE_ERROR(status) \
{cublasStatus_t result = status; \
switch (result) { \
	default: if (result != CUBLAS_STATUS_SUCCESS) {FATAL ("Other problem: " << result << " " << cudaErrorInvalidDevicePointer); throw 0;}}}

#define DEVICE_ID cuda::config::device ()

namespace utils
{
	namespace cuda
	{
		class config
		{
		public:
			config (int argc = 0, const char **argv = NULL);
	
			virtual ~config ();
		
			int device () {
				return device_id;
			}
		private:
			static int device_id;
		};
	
		template <class datatype>
		class vector 
		{
		public:
			vector (int i_n = 0, datatype *x = NULL);
	
			~vector ();
	
			datatype* ptr () {
				return data_device;
			}
		
			void resize (int i_n);
	
			int size () {
				return n;
			}
	
			void copy_to_device (int n, datatype* x);
	
			void copy_to_host (int n, datatype* x);
	
		private:
			int n;
			datatype* data_device;
		};
	} /* cuda */
} /* utils */

#endif /* end of include guard: UTILS_CUDA_CUH_AGOKALSC */
