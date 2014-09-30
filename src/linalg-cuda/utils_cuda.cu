/*!**********************************************************************
 * \file cuda.cu
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-26.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cuda_runtime.h>
#include "../extern/helper_cuda.h"
#include "utils_cuda.hpp"
#include "logger/logger.hpp"

namespace utils
{
	namespace cuda
	{
		config config_instance;
		int config::device_id;

		config::config (int argc, const char **argv) {
			device_id = findCudaDevice(argc, argv);
		}
	
		config::~config () {
		    cudaDeviceReset();
		}
	
		template <class datatype>
		vector <datatype>::vector (int i_n, datatype *x) {
			n = i_n;
			if (n != 0) {
				HANDLE_ERROR (cudaMalloc ((void**) &data_device, n * sizeof (datatype)));
				if (x) {
					copy_to_device (n, x);
				}
			}
		}
	
		template <class datatype>
		vector <datatype>::~vector () {
			if (n != 0) {
				HANDLE_ERROR (cudaFree ((datatype*) data_device));
			}
		}
	
		template <class datatype>
		void vector <datatype>::resize (int i_n) {
			if (n != 0) {
				HANDLE_ERROR (cudaFree ((datatype*) data_device));
			}
			if (i_n != 0) {
				n = i_n;
				HANDLE_ERROR (cudaMalloc ((void**) &data_device, n * sizeof (datatype)));
			}
		}
	
		template <class datatype>
		void vector <datatype>::copy_to_device (int i_n, datatype* x) {
			if (i_n != 0 && n != 0) {
				HANDLE_ERROR (cudaMemcpy(data_device, x, sizeof (datatype) * i_n, cudaMemcpyHostToDevice));
			} else if (n == 0) {
				throw 0;
			}
		}

		template <class datatype>
		void vector <datatype>::copy_to_host (int i_n, datatype* x) {
			if (i_n != 0 && n != 0) {
				HANDLE_ERROR (cudaMemcpy(x, data_device, sizeof (datatype) * i_n, cudaMemcpyDeviceToHost));
			} else if (n == 0) {
				throw 0;
			}
		}
	
		template class vector <double>;
		template class vector <float>;
		template class vector <int>;
	} /* cuda */
} /* utils */
