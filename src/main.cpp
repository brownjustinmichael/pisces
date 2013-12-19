/*!***********************************************************************
 * \file main.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "bases/messenger.hpp"
#include "one_d/element_one_d.hpp"
#include "two_d/element_two_d.hpp"
#include "config.hpp"

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
 * \namespace one_d::chebyshev
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
	
	// Initialize messenger
	bases::messenger process_messenger (&argc, &argv, 2);

	id = process_messenger.get_id ();
	n_elements = process_messenger.get_np ();

	log_config log_config_instance (&argc, &argv, id);
	
	io::parameters <double> params ("../input/parameters.txt");

	int n = params.gridpoints / n_elements + 1;
	double position_m0 = -1.0 + 2.0 / n_elements * id;
	double position_mm = -1.0 + 2.0 / n_elements * (id + 1);
	double position_n0 = -1.0;
	double position_nn = 1.0;
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
	
	int m = params.gridpoints;
	
	bases::axis horizontal_axis (m, 0, position_n0, 0, position_nn);
	bases::axis vertical_axis (n, excess_0, position_m0, excess_n, position_mm);
		
	// one_d::chebyshev::advection_diffusion_element <double> element (&vertical_axis, name, params, &process_messenger, 0x00);
	two_d::fourier::chebyshev::advection_diffusion_element <double> element (&horizontal_axis, &vertical_axis, name, params, &process_messenger, 0x00);
	
	try {
		element.run ();
	} catch (...) {
		FATAL ("Fatal error occurred. Check log.");
		return 1;
	}
	
	INFO ("Main complete.");
	
	return 0;
}
