/*!***********************************************************************
 * \file main.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "config.hpp"
#include "one_d/element_one_d.hpp"
#include "one_d/boundary_one_d.hpp"

int config::n_loggers = 0;
int config::n_appenders = 0;
int config::severity = 4;  // The default logging severity is 4, errors and fatal messages only.
std::string config::config_file = "../input/Log4cxxConfig.xml";
int n = 256; //!< The number of grid points
double scale = 1.0;
double sigma = 0.1;

std::vector<double> initial_conditions;
std::vector<double> initial_position;

#ifdef __APPLE__

#include <log4cxx/patternlayout.h>

// Set up the logs

log4cxx::LayoutPtr config::console_layout = new log4cxx::PatternLayout ("%-5p %c{2} (%C::%M %L) - %m%n");
log4cxx::LayoutPtr config::layout = new log4cxx::PatternLayout ("%d %-5p [%t:%x] %C{2} (%C::%M %L) - %m%n");
std::vector<log4cxx::LoggerPtr> config::loggers;
std::vector<log4cxx::AppenderPtr> config::appenders;

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
		
	config::make_main ();
	
	MTRACE ("Beginning main...");
	
	initial_position.resize (n + 1);
	initial_conditions.resize (n + 1);
	
	double pioN = std::acos (-1.0) / (n / 2 - 1);
	for (int i = 0; i < n / 2; ++i) {
		initial_position [i] = std::cos (pioN * i) + 1.0;
		initial_position [i + n / 2] = std::cos (pioN * i) - 1.0;
		initial_conditions [i] = scale * std::exp (- (initial_position [i] - 0.0) * (initial_position [i] - 0.0) / 2.0 / sigma / sigma) - scale * std::exp (- 1.0 / 2.0 / sigma / sigma);
		initial_conditions [i + n / 2] = scale * std::exp (- (initial_position [i + n / 2] - 0.0) * (initial_position [i + n / 2] - 0.0) / 2.0 / sigma / sigma) - scale * std::exp (- 1.0 / 2.0 / sigma / sigma);
	}
	initial_position [n] = -2.0;
	initial_conditions [n] = 0.0;
	
	one_d::advection_diffusion_element element_1 ("_1_", 0.0, 0.0, n / 2, &initial_position [0], &initial_conditions [0], 0x00);
	one_d::advection_diffusion_element element_2 ("_2_", 0.0, 0.0, n / 2, &initial_position [n / 2], &initial_conditions [n / 2], 0x00);
	element_1.add_boundary (one_d::boundary::make_unique (0.0, element_1 (rhs)));
	element_1.add_boundary (one_d::boundary::make_unique (0.0, element_1 (rhs, n / 2 - 1)));
	element_2.add_boundary (one_d::boundary::make_unique (0.0, element_2 (rhs)));
	element_2.add_boundary (one_d::boundary::make_unique (0.0, element_2 (rhs, n / 2 - 1)));
	// element_1.add_boundary (one_d::boundary::make_unique (0.5, element_1 (rhs, n / 2 - 1), 0.5, element_2 (rhs)));

	MTRACE ("main: Entering main loop.");
	
	for (i = 0; i < 500; ++i) {
		MTRACE ("main: Beginning timestep...");
		MINFO ("main: Timestep: " << i);

		element_1.calculate ();
		element_2.calculate ();
		element_1.execute_boundaries ();
		element_2.execute_boundaries ();
		element_1.update ();
		element_2.update ();
		
		MTRACE ("main: Timestep " << i << " complete.");
	}

	MTRACE ("main: End of main.");

	return 0;
}