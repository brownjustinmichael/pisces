// 
//! \file main.cpp
//  spectral element
//  
//  Created by Justin Brown on 2013-03-22.
//  Copyright 2013 Justin Brown. All rights reserved.
// 

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <cmath>
#include <fftw3.h>
#include "config.hpp"
#include "diffusion/diffusion.hpp"
#include "boundary/boundary.hpp"
#include "io/io.hpp"
#include "io/exceptions.hpp"

#define N 16

log4cxx::LoggerPtr config::logger (log4cxx::Logger::getLogger ("main"));
log4cxx::LevelPtr config::levels [6];

/*! \mainpage
 *
 *  \author Justin Brown
 *  \author Ryan Moll
 *
 *  \section Introduction
 *
 *  The goal of this project is to set up a code designed to do 2D Anelastic simulations using a spectral element scheme. Currently, the aim is do get the solution to the 1D advection-diffusion equation for a constant background flow.
 *
 *  Possible further reaching goals: 3D, pseudo-incompressible
 */

int main (int argc, char const *argv[])
{	
	int i, j;
	
	// Here, we're setting up the log. The levels array allows us to easily convert from integer to logging severity, which is specified by the -D flag at the command line, e.g. -D3 specifies that only warnings, errors, and fatal messages will be logged.
	log4cxx::xml::DOMConfigurator::configure("../input/Log4cxxConfig.xml");
	
	config::levels [0] = log4cxx::Level::getTrace ();
	config::levels [1] = log4cxx::Level::getDebug ();
	config::levels [2] = log4cxx::Level::getInfo ();
	config::levels [3] = log4cxx::Level::getWarn ();
	config::levels [4] = log4cxx::Level::getError ();
	config::levels [5] = log4cxx::Level::getFatal ();

	config::logger->setLevel (config::levels [4]); // The default logging severity is 4, errors and fatal messages only.
	
	// The program runs through the execution flags.
	while ((argc > 1) && (argv [1] [0] == '-')) {
		switch (argv [1] [1]) {
			// Debug switch
			case 'D':
				config::logger->setLevel (config::levels [atoi (&(argv [1] [2]))]);
				break;
		}
		--argc;
		++argv;
	}
	
	LOG4CXX_TRACE (config::logger, "Beginning main...");
	
	double timestep = 0.001;
	
	std::vector<double> velocity (N, 0.0);
	std::vector<double> position (N, 0.0);
	
	double pioN = std::acos (-1.0) / N;
	
	for (i = 0; i < N; ++i) {
		position [i] = std::cos (pioN * i);
	}
	
	std::vector<double *> data_ptrs (2);
	data_ptrs [0] = &position [0];
	data_ptrs [1] = &velocity [0];
	
	velocity [0] = 2.0;
	velocity [1] = 0.0;
	velocity [2] = -1.0;
	velocity [3] = 0.0;
	
	io::incremental_output_stream_1D cheb_stream ("../output/test_cheb", ".dat", 4, new io::header, N, 1, &data_ptrs [1]);
	io::incremental_output_stream_1D angle_stream ("../output/test_angle", ".dat", 4, new io::header, N, 2, &data_ptrs [0]);
	io::simple_output_stream_1D failsafe_dump ("_dump.dat", N, 2, &data_ptrs [0]);
		
	diffusion::cheb_1D diffusion_plan (1., N, &velocity [0]);
	
	boundary::fixed_cart_1D upper_bound (&velocity [0], 0.0);
	boundary::fixed_cart_1D lower_bound (&velocity [N - 1], 0.0);
	
	fftw_plan fourier_plan = fftw_plan_r2r_1d (N, &velocity [0], &velocity [0], FFTW_REDFT00, FFTW_ESTIMATE);

	LOG4CXX_TRACE (config::logger, "main: Entering main loop.");

	// Output in Chebyshev space
	try {
		cheb_stream.output ();
	} catch (io::exceptions::file_exception &io_exception) {
		LOG4CXX_ERROR (config::logger, "Unable to print to file, outputting failsafe dump to _dump.dat");
		failsafe_dump.output ();
		exit (EXIT_FAILURE);
	}
	// Transform forward
	fftw_execute (fourier_plan);
	// Output in angle space
	angle_stream.output ();
	// Transform backward
	fftw_execute (fourier_plan);
	// We rescale the values to account for the factor of 2(N-1) that occurs during FFT
	for (j = 0; j < N; ++j) {
		velocity [j] /= (2 * (N - 1));
	}
	
	for (i = 0; i < 10; ++i) {
		LOG4CXX_TRACE (config::logger, "main: Beginning timestep...");
		LOG4CXX_INFO (config::logger, "main: Timestep: " << i);

		// Output in Chebyshev space
		try {
			cheb_stream.output ();
		} catch (io::exceptions::file_exception &io_exception) {
			LOG4CXX_ERROR (config::logger, "Unable to print to file, outputting failsafe dump to _dump.dat");
			failsafe_dump.output ();
			exit (EXIT_FAILURE);
		}
		// Calculate the diffusion in Chebyshev space
		diffusion_plan.execute (timestep);
		
		// Transform forward
		fftw_execute (fourier_plan);
		
		// Apply boundary conditions
		upper_bound.execute (timestep);
		lower_bound.execute (timestep);

		// Output in angle space
		angle_stream.output ();
		
		// Transform backward
		fftw_execute (fourier_plan);
		
		// We rescale the values to account for the factor of 2(N-1) that occurs during FFT
		for (j = 0; j < N; ++j) {
			velocity [j] /= (2 * (N - 1));
		}

		LOG4CXX_TRACE (config::logger, "main: Timestep " << i << " complete.");
	}

	LOG4CXX_TRACE (config::logger, "main: End of main.");

	return 0;
}