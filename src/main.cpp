/*!***********************************************************************
 * \file main.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "config.hpp"
#include "boundary/boundary.hpp"
#include "element/element.hpp"

#define N 128

int config::severity = 4;  // The default logging severity is 4, errors and fatal messages only.

#ifdef __APPLE__

// Set up the log

log4cxx::LoggerPtr config::logger (log4cxx::Logger::getLogger ("main"));

#endif // __APPLE__

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

int main (int argc, char const *argv[])
{	
	int i;
	
	// The program runs through the execution flags.
	while ((argc > 1) && (argv [1] [0] == '-')) {
		switch (argv [1] [1]) {
			// Debug switch
			case 'D':
				config::severity = atoi (&(argv [1] [2])); break;
		}
		--argc;
		++argv;
	}
	
#ifdef __APPLE__
	
	log4cxx::xml::DOMConfigurator::configure("../input/Log4cxxConfig.xml");
	config::logger->setLevel (config::int_to_severity (config::severity));
	
#endif // __APPLE__
	
	TRACE ("Beginning main...");
	
	element::diffusion_element_1D main_element (N, boundary::fixed_upper | boundary::fixed_lower);

	TRACE ("main: Entering main loop.");
	
	for (i = 0; i < 10; ++i) {
		TRACE ("main: Beginning timestep...");
		INFO ("main: Timestep: " << i);

		main_element.update ();

		TRACE ("main: Timestep " << i << " complete.");
	}

	TRACE ("main: End of main.");

	return 0;
}