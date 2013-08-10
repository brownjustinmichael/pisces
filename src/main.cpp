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
#include "one_d/master_one_d.hpp"

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
	
	// Initialize CUBLAS
#ifdef CUBLAS
	cublasInit ();
#endif

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
	
	const int n_total = inputParams ["n_elements"].asInt;
	int n;
	int n_masters;
	int index;
	std::vector <int> n_grid (n_total, inputParams ["gridpoints"].asInt / n_total);
	std::vector <double> position_grid (n_total + 1, 0.0);
	position_grid [0] = -1.0;
	for (int i = 1; i < n_total + 1; ++i) {
		position_grid [i] = position_grid [i - 1] + 2.0 / n_total;
	}
	std::vector <std::string> name_grid (n_total);
	for (int i = 0; i < n_total; ++i) {
		name_grid [i] = std::to_string (i + 1);
	}
	
	if (process_messenger.get_id () >= n_total % process_messenger.get_np ()) {
		n = n_total / process_messenger.get_np ();
		index = n * process_messenger.get_id () + n_total % process_messenger.get_np ();
	} else {
		n = n_total / process_messenger.get_np () + 1;
		index = n * process_messenger.get_id ();
	}

	if (process_messenger.get_id () >= n_total) {
		return 0;
	} else {
		if (process_messenger.get_np () > n_total) {
			n_masters = n_total;
		} else {
			n_masters = process_messenger.get_np ();
		}
	}
	
	one_d::master <one_d::chebyshev::advection_diffusion_element> master_process (process_messenger.get_id (), n_masters, "../input/parameters.txt", n, &n_grid [index], &position_grid [index], &process_messenger);
	
	if (process_messenger.get_id () != 0) {
		MTRACE ("Adding boundary to " << 0 << " at 0 at processor " << process_messenger.get_id () - 1);
		master_process.add_boundary (0, one_d::edge_0, 1, 2, process_messenger.get_id () - 1);
	}
	if (process_messenger.get_id () != n_masters - 1) {
		MTRACE ("Adding boundary to " << n - 1 << " at n - 1 at processor " << process_messenger.get_id () + 1);
		master_process.add_boundary (n - 1, one_d::edge_n, 2, 1, process_messenger.get_id () + 1);
	}
	
	master_process.run ();
	
	return 0;
}
