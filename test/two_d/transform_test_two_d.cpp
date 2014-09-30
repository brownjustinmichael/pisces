/*!**********************************************************************
 * \file transform_test.cpp
 * /Developer/NVIDIA/CUDA-5.5/samples/7_CUDALibraries/simpleCUFFT
 * 
 * Created by Justin Brown on 2014-03-01.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "messenger/messenger.hpp"
#include "io/io.hpp"
#include "element/element_one_d.hpp"
#include "logger/logger.hpp"
#include <fftw3.h>

#define TINY 1.e-4

using namespace two_d::fourier::cosine;

int main (int argc, char *argv[])
{
	log_config::configure (&argc, &argv);

	int n = 100;
	std::vector <double> init (n * n), init_copy (n * n), fftw_copy (n * n + 1);
	
	// {
	// 	scalar = 1.0 / std::sqrt (n);
	// 					
	// 	iodim.n = n;
	// 	major_iodim.is = 1;
	// 	major_iodim.os = 1;
	// 	
	// 	plans_float.resize (threads);
	// 	
	// 	if (!(flags & inverse)) {
	// 		iodim.is = m;
	// 		iodim.os = 2 * m;
	// 		int index = 0;
	// 		for (int i = 0; i < threads; ++i) {
	// 			major_iodim.n = m / threads + (i < (m % threads)? 1: 0);
	// 			plans [i] = fftw_plan_guru_split_dft_r2c (1, &iodim, 1, &major_iodim, data_in + index, data_out + index, data_out + m + index, FFTW_MEASURE);
	// 			index += major_iodim.n;
	// 		}
	// 	} else {
	// 		iodim.is = 2 * m;
	// 		iodim.os = m;
	// 		int index = 0;
	// 		for (int i = 0; i < threads; ++i) {
	// 			major_iodim.n = m / threads + (i < (m % threads)? 1: 0);
	// 			plans [i] = fftw_plan_guru_split_dft_c2r (1, &iodim, 1, &major_iodim, data_in + index, data_in + m + index, data_out + index, FFTW_MEASURE);
	// 			index += major_iodim.n;
	// 		}
	// 	}
	// }
	
	// {
	// 	scalar = 1.0;
	// 	if (m > 1 && !(flags & ignore_m)) {
	// 		scalar /= std::sqrt (2.0 * (m - 1));
	// 	}
	// 	
	// 	fftw_r2r_kind kind = FFTW_REDFT00;
	// 
	// 	plans.resize (threads);
	// 	
	// 	int index = 0;
	// 	for (int i = 0; i < threads; ++i) {
	// 		int nn = (2 * (n / 2 + 1)) / threads + (i < (2 * (n / 2 + 1) % threads)? 1: 0);
	// 		plans [i] = fftw_plan_many_r2r (1, &m, nn, data_in + index * m, NULL, 1, m, data_out + index * m, NULL, 1, m, &kind, FFTW_MEASURE);
	// 		index += nn;
	// 	}
	// }
	fftw_plan fourier_plan = fftw_plan_r2r_1d (n, &fftw_copy [0], &fftw_copy [0], FFTW_REDFT00, FFTW_ESTIMATE);
	
	utils::messenger process_messenger (&argc, &argv);
	io::parameters config;
	
	config ["time.max"] = 0.0;
	config ["time.cfl"] = 0.0;
	config ["time.alpha"] = 0.0;
	config ["grid.z.width"] = 0.0;
	config ["init.scale"] = 0.0;
	config ["init.mean"] = 0.0;
	config ["init.sigma"] = 0.0;
	config ["velocity.advection"] = 0.0;
	config ["velocity.diffusion"] = 0.0;
	config ["output.file"] = "";
	config ["output.every"] = 0;
	config ["parallel.transform.threads"] = 0;
	
	bases::axis vertical_axis (n, -1.0, 1.0);
	bases::axis horizontal_axis (m, -1.0, 1.0);
	
	boussinesq_element <double> element (&horizontal_axis, &vertical_axis, 0, config, &process_messenger, 0x00);
	
	srand (1);
	for (int i = 0; i < n; ++i) {
		init [i] = rand () % 100;
		init_copy [i] = init [i];
		fftw_copy [i] = init [i];
	}
		
	element.initialize (0, &init [0]);
	
	element.transform (forward_vertical);
	
	fftw_execute (fourier_plan);
	
	for (int i = 0; i < n; ++i) {
		if (abs (element.ptr (0) [i] - fftw_copy [i] * 1.0 / std::sqrt (2.0 * (n - 1))) > TINY) {
			FATAL ("Forward transform discrepancy in 1D Chebyshev Advection Diffusion Element");
			return 1;
		}
	}
	
	element.transform (inverse_vertical);
	
	for (int i = 0; i < n; ++i) {
		if (abs (element.ptr (0) [i] - init_copy [i]) > TINY) {
			FATAL ("Inverse transform discrepancy in 1D Chebyshev Advection Diffusion Element");
			return 1;
		}
	}
	return 0;
}