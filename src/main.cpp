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
#include <memory>
#include <map>
#include <vector>
#include "config.hpp"
#include "one_d/element_one_d.hpp"
#include "one_d/boundary_one_d.hpp"
#include "one_d/master_one_d.hpp"

// #include "mpi.h"

int config::n_loggers = 0;
int config::n_appenders = 0;
int config::severity = 4;  // The default logging severity is 4, errors and fatal messages only.
std::string config::config_file = "../input/Log4cxxConfig.xml";
unsigned int n = 66; //!< The number of grid points
double scale = 1.0;
double sigma = 0.1;

std::vector<double> initial_conditions;
std::vector<double> initial_position;

#ifdef __APPLE__

#include <log4cxx/patternlayout.h>

// Set up the logs

log4cxx::LayoutPtr config::console_layout = new log4cxx::PatternLayout ("%-5p %c{2}: %C (%M %L) - %m%n");
log4cxx::LayoutPtr config::layout = new log4cxx::PatternLayout ("%d %-5p %c{2}: %C (%M %L) - %m%n");
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
	// int id;
	// int p;
	// 
	// // Initialize mpi
	// MPI::Init (argc, argv);
	// 
	// // Get the number of processes
	// p = MPI::COMM_WORLD.Get_size ();
	
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
	
	std::vector <int> n_grid (2, 32);
	std::vector <double> position_grid (3, 0.0);
	position_grid [0] = 2.0;
	position_grid [2] = -2.0;
	std::vector <std::string> name_grid (2);
	name_grid [0] = "1";
	name_grid [1] = "2";
	
	one_d::master <one_d::chebyshev::advection_diffusion_element, one_d::diffusive_boundary> master_process (1, "../input/parameters.txt", 2, &n_grid [0], &position_grid [0], &name_grid [0]);
	
	master_process.run ();
	
	// if (id == 0) {
		// Read experiment parameters out of text file
		// std::string filename = "../input/parameters.txt";
		// io::parameter_map inputParams;
		// io::read_params_txt parameters (filename);
		// inputParams = parameters.load_params();
		// 
		// // Inialize some experiment variables
		// const int tsteps = inputParams["timesteps"].asInt;
		// 	
		// MTRACE ("Beginning main..." << n);
	// }
	
	// one_d::chebyshev::advection_diffusion_element element_1 (n / 2, 0.0, 2.0, "1", inputParams, 0x00);
	// one_d::chebyshev::advection_diffusion_element element_2 (n / 2, -2.0, 0.0, "2", inputParams, 0x00);
	// 
	// std::vector <double> buffer_1 (n / 2), buffer_2 (n / 2);
	// 
	// element_1.add_boundary (std::make_shared <one_d::diffusive_boundary> (one_d::diffusive_boundary (linked_n, inputParams["diffusion_coeff"].asDouble, position, velocity, rhs, &buffer_1 [0], &buffer_2 [0])));
	// element_2.add_boundary (std::make_shared <one_d::diffusive_boundary> (one_d::diffusive_boundary (linked_0, inputParams["diffusion_coeff"].asDouble, position, velocity, rhs, &buffer_2 [0], &buffer_1 [0])));
	// 
	// MTRACE ("main: Entering main loop.");
	// 
	// for (int i = 0; i < tsteps; ++i) {
	// 	MTRACE ("main: Beginning timestep...");
	// 	MINFO ("main: Timestep: " << i);
	// 	
	// 	try {
	// 		element_1.calculate ();
	// 		element_2.calculate ();
	// 		element_1.send ();
	// 		element_2.send ();
	// 		element_1.recv ();
	// 		element_2.recv ();
	// 		element_1.execute_boundaries ();
	// 		element_2.execute_boundaries ();
	// 		element_1.update ();
	// 		element_2.update ();
	// 	} catch (...) {
	// 		element_1.failsafe ();
	// 		element_2.failsafe ();
	// 		exit (EXIT_FAILURE);
	// 	}
	// 	
	// 	MTRACE ("main: Timestep " << i << " complete.");
	// }
	// 
	// MTRACE ("main: End of main.");
	
	// MPI::Finalize ();

	return 0;
}