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
#include "config.hpp"
#include "diffusion/diffusion.hpp"
#include "io/io.hpp"

#define N 8

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
	int i;
	
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
					
	std::vector<double> velocity (N, 0.0);
	std::vector<double *> data_ptrs (2);
	data_ptrs [0] = &velocity[0];
	data_ptrs [1] = &velocity[0];
	
	velocity [0] = 1.0;
	velocity [1] = 1.0;
	velocity [2] = 1.0;
	velocity [3] = 1.0;
	
	io::incremental_output_stream_1D test_stream ("../data/test", ".dat", 4, new io::header, N, 1, &data_ptrs [0]);
		
	diffusion::cheb_1D diffusion_operation (N);
	diffusion::plan diffusion_plan (&diffusion_operation, 2, &data_ptrs [0]);

	LOG4CXX_TRACE (config::logger, "main: Entering main loop.");

	for (i = 0; i < 10; ++i) {
		LOG4CXX_TRACE (config::logger, "main: Beginning timestep...");
		LOG4CXX_INFO (config::logger, "main: Timestep: " << i);

		test_stream.output ();
		diffusion_plan.execute (0.01);

		LOG4CXX_TRACE (config::logger, "main: Timestep " << i << " complete.");
	}

	LOG4CXX_TRACE (config::logger, "main: End of main.");

	return 0;
}