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

log4cxx::LoggerPtr config::logger (log4cxx::Logger::getLogger ("main"));
log4cxx::LevelPtr config::levels [6];

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
	
	// Here, we're setting up the log. The levels array allows us to easily convert from integer to logging severity, which is specified by the -D flag at the command line, e.g. -D3 specifies that only warnings, errors, and fatal messages will be logged.
	log4cxx::xml::DOMConfigurator::configure("../input/Log4cxxConfig.xml");
	config::logger->setLevel (config::levels [4]); // The default logging severity is 4, errors and fatal messages only.
	
	// The program runs through the execution flags.
	while ((argc > 1) && (argv [1] [0] == '-')) {
		switch (argv [1] [1]) {
			// Debug switch
			case 'D':
				config::logger->setLevel (config::int_to_severity (atoi (&(argv [1] [2]))));
				break;
		}
		--argc;
		++argv;
	}
	
	TRACE ("Beginning main...");
	
	element::diffusion_element main_element (N, boundary::fixed_upper | boundary::fixed_lower);

	TRACE ("main: Entering main loop.");
	
	for (i = 0; i < 100; ++i) {
		TRACE ("main: Beginning timestep...");
		INFO ("main: Timestep: " << i);

		main_element.update ();

		TRACE ("main: Timestep " << i << " complete.");
	}

	TRACE ("main: End of main.");

	return 0;
}