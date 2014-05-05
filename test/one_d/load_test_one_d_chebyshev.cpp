/*!**********************************************************************
 * \file load_test.cpp
 * /Developer/NVIDIA/CUDA-5.5/samples/7_CUDALibraries/simpleCUFFT
 * 
 * Created by Justin Brown on 2014-03-01.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "../../src/bases/messenger.hpp"
#include "../../src/utils/io.hpp"
#include "../../src/one_d/element_one_d.hpp"

using namespace one_d::chebyshev;

int main (int argc, char const *argv[])
{
	bases::messenger process_messenger (&argc, &argv, 2);
	io::parameters config ("../input/config.yaml");
	bases::axis horizontal_axis (n, position_n0, position_nn);

	advection_diffusion_element <double> element (&vertical_axis, name, config, &process_messenger, 0x00);	

	return 0;
}