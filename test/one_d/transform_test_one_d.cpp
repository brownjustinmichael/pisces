/*!**********************************************************************
 * \file transform_test.cpp
 * /Developer/NVIDIA/CUDA-5.5/samples/7_CUDALibraries/simpleCUFFT
 * 
 * Created by Justin Brown on 2014-03-01.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "../../src/bases/messenger.hpp"
#include "../../src/utils/io.hpp"
#include "../../src/one_d/element_one_d.hpp"
#include "../../src/config.hpp"
#include <fftw3.h>

#define TINY 1.e-4

using namespace one_d::chebyshev;

int main (int argc, char *argv[])
{
	log_config::configure (&argc, &argv);

	int n = 100;
	std::vector <double> init (n), init_copy (n), fftw_copy (n + 1);
	fftw_plan fourier_plan = fftw_plan_r2r_1d (n, &fftw_copy [0], &fftw_copy [0], FFTW_REDFT00, FFTW_ESTIMATE);
	
	bases::messenger process_messenger (&argc, &argv);
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
	
	advection_diffusion_element <double> element (&vertical_axis, 0, config, &process_messenger, 0x00);
	
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
	
	INFO ("One_d transform test passed.");
	return 0;
}