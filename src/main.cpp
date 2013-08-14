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

#include "mpi.h"

int config::n_loggers = 0;
int config::n_appenders = 0;
int config::severity = 4;  // The default logging severity is 4, errors and fatal messages only.
std::string config::config_file = "../input/Log4cxxConfig.xml";

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
int main (int argc, char *argv[])
{
	// Initialize messenger
	utils::messenger process_messenger (&argc, &argv);

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
	
	config::make_main (process_messenger.get_id ());
	
	io::parameter_map inputParams;
	io::read_params_txt parameters ("../input/parameters.txt");
	inputParams = parameters.load_params();
	
	int id = process_messenger.get_id ();
	int n_elements = process_messenger.get_np ();
	int n = inputParams ["gridpoints"].asInt / n_elements;
	double position_0 = -1.0 + 2.0 / n_elements * id;
	double position_n = -1.0 + 2.0 / n_elements * (id + 1);
	int excess_0;
	int excess_n;
	if (id == 0) {
		excess_0 = 0;
	} else {
		excess_0 = 1;
	}
	if (id == n_elements - 1) {
		excess_n = 0;
	} else {
		excess_n = 1;
	}
	int name = id;
	
	one_d::chebyshev::advection_diffusion_element element (n, position_0, position_n, excess_0, excess_n, 0, name, inputParams, &process_messenger, 0x00);
	
	if (id != 0) {
		MTRACE ("Adding boundary to " << name << " at 0 at processor " << id - 1);
		element.add_boundary (one_d::edge_0, 1, 2, id - 1, 0);
	}
	if (id != n_elements - 1) {
		MTRACE ("Adding boundary to " << name << " at n - 1 at processor " << id + 1);
		element.add_boundary (one_d::edge_n, 2, 1, id + 1, 0);
	}
	
	double t_timestep;
	for (int i = 0; i < inputParams ["timesteps"].asInt; ++i) {
		MINFO ("Timestep " << i);
		element.calculate ();
		element.output ();
		element.execute_boundaries ();
		t_timestep = element.calculate_timestep ();
		process_messenger.min (&t_timestep);
		MTRACE ("Updating...");
		for (int k = 0; k < 2; ++k) {
			element.send_positions ();
			element.recv_positions ();
			element.attempt_update ();
			element.calculate_bounds ();
			element.send_bounds ();
			element.recv_bounds ();
		}
		element.attempt_update ();
		element.update ();
		element.update_timestep (t_timestep);
	}
	
	return 0;
}
