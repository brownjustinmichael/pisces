/*!**********************************************************************
 * \file main_cuda.cpp
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-14.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "utils/utils_cublas.hcu"
#include "config.hpp"

/*!*******************************************************************
 * \mainpage
 *
 * \author Justin Brown
 * \author Ryan Moll
 *
 * \section Introduction
 *
 * The goal of this project is to set up a code designed to do 2D 
 * Anelastic simulations using a spectral element scheme. Currently, 
 * the aim is do get the solution to the 1D advection-diffusion 
 * equation for a constant background flow.
 *
 * Possible further reaching goals: 3D, pseudo-incompressible
 *********************************************************************/

/*!*******************************************************************
 * \namespace bases
 * 
 * \brief A namespace containing the base classes of the code
 *********************************************************************/

/*!*******************************************************************
 * \namespace io
 * 
 * \brief A namespace containing all the input and output classes of the code
 *********************************************************************/

/*!*******************************************************************
 * \namespace utils
 * 
 * \brief A namespace containing the various utilities needed in the code, such as linear algebra
 *********************************************************************/

/*!*******************************************************************
 * \namespace one_d
 * 
 * \brief A namespace containing all the 1D pieces of the code
 *********************************************************************/

/*!*******************************************************************
 * \namespace chebyshev
 * 
 * \brief A namespace containing the 1D Chebyshev pieces of the code
 *********************************************************************/

/*!*******************************************************************
 * \namespace two_d
 * 
 * \brief A namespace containing all the 2D pieces of the code
 *********************************************************************/

/*!*******************************************************************
 * \brief The main call
 * 
 * \param argc The integer number of command line arguments
 * \param argv The character array of command line arguments
 *********************************************************************/

void test () {
	// io::parameter_map inputParams;
	// io::read_params_txt parameters ("../input/parameters.txt");
	// inputParams = parameters.load_params();
	// 
	// 
	// int n = inputParams ["gridpoints"].asInt / n_elements;
	// double position_0 = -1.0 + 2.0 / n_elements * id;
	// double position_n = -1.0 + 2.0 / n_elements * (id + 1);
	// int excess_0;
	// int excess_n;
	// if (id == 0) {
	// 	excess_0 = 0;
	// } else {
	// 	excess_0 = 1;
	// }
	// if (id == n_elements - 1) {
	// 	excess_n = 0;
	// } else {
	// 	excess_n = 1;
	// }
	// int name = id;
	
	float a [100];
		
	for (int i = 0; i < 100; ++i) {
		a [i] = i;
	}
	
	utils::cuda::vector <float> dev_a (100, a), dev_b (100);
	
	utils::cuda::scale (100, (float) 2.0, &dev_a);
	
	dev_a.copy_to_host (100, a);
	
	for (int i = 0; i < 100; ++i) {
		printf ("Final %d: %f\n", i, a [i]);
	}
	
	// one_d::chebyshev::cuda_element element (n, position_0, position_n, excess_0, excess_n, name, inputParams, &process_messenger, 0x00);
	// 
	// if (id != 0) {
	// 	TRACE ("Adding boundary to " << name << " at 0 at processor " << id - 1);
	// 	element.add_boundary (one_d::edge_0, 1, 2, id - 1);
	// }
	// if (id != n_elements - 1) {
	// 	TRACE ("Adding boundary to " << name << " at n - 1 at processor " << id + 1);
	// 	element.add_boundary (one_d::edge_n, 2, 1, id + 1);
	// }
	// 
	// element.send_positions ();
	// element.recv_positions ();
	// 
	// element.run ();
	// 
	// INFO ("Main complete.");
	// 
}