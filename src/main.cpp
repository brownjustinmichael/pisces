/*!***********************************************************************
 * \file main.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cmath>
#include <vector>
#include "config.hpp"
#include "one_d/element_one_d.hpp"
#include "one_d/boundary_one_d.hpp"

int config::severity = 4;  // The default logging severity is 4, errors and fatal messages only.
int n = 256; //!< The number of grid points
double scale = 1.0;
double sigma = 0.1;

std::vector<double> initial_conditions;
std::vector<double> initial_position;

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

/*!*******************************************************************
 * \brief The main call
 * 
 * \param argc The integer number of command line arguments
 * \param argv The character array of command line arguments
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
	
	initial_position.resize (n);
	initial_conditions.resize (n);
	
	double pioN = std::acos (-1.0) / n;
	for (int i = 0; i < n; ++i) {
		initial_position [i] = std::cos (pioN * i);
		initial_conditions [i] = scale * std::exp (- initial_position [i] * initial_position [i] / 2.0 / sigma / sigma) - scale * std::exp (- 1.0 / 2.0 / sigma / sigma);
	}
	
	one_d::advection_diffusion_element element_1 (n, &initial_conditions [0], 0x00);
	// one_d::advection_diffusion_element element_2 (n / 2, &intial_conditions [n / 2], 0x00);
	// add_boundary (boundary::make_unique (0.0, element_1 [velocity], 0.0, &(scalars [velocity] [n - 1])));
	// add_boundary (boundary::make_unique (0.0, &(scalars [rhs] [0]), 0.0, &(scalars [rhs] [n - 1])));

	TRACE ("main: Entering main loop.");
	
	for (i = 0; i < 500; ++i) {
		TRACE ("main: Beginning timestep...");
		INFO ("main: Timestep: " << i);

		element_1.calculate ();
		element_1.execute_boundaries ();
		element_1.update ();
		
		TRACE ("main: Timestep " << i << " complete.");
	}

	TRACE ("main: End of main.");

	return 0;
}