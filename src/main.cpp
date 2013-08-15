/*!***********************************************************************
 * \file main.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "config.hpp"
#include "one_d/element_one_d.hpp"

#ifdef _MPI
#include "mpi.h"
#endif // _MPI

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
 * \brief A namespace containing all the input and output classes of 
 * the code
 *********************************************************************/

/*!*******************************************************************
 * \namespace utils
 * 
 * \brief A namespace containing the various utilities needed in the 
 * code, such as linear algebra
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
	int id = 0, n_elements = 1;
	
#if _MPI
	// Initialize messenger
	utils::messenger process_messenger (&argc, &argv);

	id = process_messenger.get_id ();
	n_elements = process_messenger.get_np ();
#endif // _MPI

	log_config::update_name (id);

	// The program runs through the execution flags.
	while ((argc > 1) && (argv [1] [0] == '-')) {
		switch (argv [1] [1]) {
			// Debug switch
			case 'D':
				log_config::update_severity (atoi (&(argv [1] [2])));
				break;
		}
		--argc;
		++argv;
	}
	
	io::parameter_map inputParams;
	io::read_params_txt parameters ("../input/parameters.txt");
	inputParams = parameters.load_params();
	

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
	
	one_d::chebyshev::advection_diffusion_element element (n, position_0, position_n, excess_0, excess_n, name, inputParams, &process_messenger, 0x00);
	
	if (id != 0) {
		TRACE ("Adding boundary to " << name << " at 0 at processor " << id - 1);
		element.add_boundary (one_d::edge_0, 1, 2, id - 1);
	}
	if (id != n_elements - 1) {
		TRACE ("Adding boundary to " << name << " at n - 1 at processor " << id + 1);
		element.add_boundary (one_d::edge_n, 2, 1, id + 1);
	}

	element.send_positions ();
	element.recv_positions ();
	
	element.run ();
	
	INFO ("Main complete.");
	
	return 0;
}
