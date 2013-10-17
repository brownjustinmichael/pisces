/*!***********************************************************************
 * \file main.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "config.hpp"
#include "two_d/element_two_d.hpp"
#include "two_d/transform_two_d.hpp"

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
	

	int n = inputParams ["gridpoints"].asInt / n_elements + 1;
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

	// if (id != 0) {
	// 	TRACE ("Adding boundary to " << name << " at 0 at processor " << id - 1);
	// 	process_messenger.add_boundary (one_d::edge_0, id - 1);
	// }
	// if (id != n_elements - 1) {
	// 	TRACE ("Adding boundary to " << name << " at n - 1 at processor " << id + 1);
	// 	process_messenger.add_boundary (one_d::edge_n, id + 1);
	// }
	
	int m = n;
	
	bases::axis horizontal_axis (n, excess_0, position_0, excess_n, position_n);
	bases::axis vertical_axis (m, excess_0, position_0, excess_n, position_n);
	
	// two_d::fourier::chebyshev::advection_diffusion_element <double> element (&horizontal_axis, &vertical_axis, name, inputParams, &process_messenger, 0x00);
	
	bases::fourier::grid <double> horizontal_grid (&horizontal_axis);
	bases::chebyshev::grid <double> vertical_grid (&vertical_axis);
	
	std::vector <int> cell_n (n * m);
	std::vector <int> cell_m (n * m);
	std::vector <double> position_h (n * m);
	std::vector <double> position_m (n * m);
	std::vector <double> data (m * (n + 2));
	
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			position_h [i * m + j] = horizontal_grid [i];
			position_m [i * m + j] = vertical_grid [j];
			cell_n [i * m + j] = i;
			cell_m [i * m + j] = j;
			data [i * m + j] = exp (- (vertical_grid [j] * vertical_grid [j] + horizontal_grid [i] * horizontal_grid [i]) * 16.0);
		}
	}
	
	io::simple_output <double> out ("output.dat", n * m);
	out.append (&cell_n [0]);
	out.append (&cell_m [0]);
	out.append (&position_h [0]);
	out.append (&position_m [0]);
	out.append (&data [0]);
		
	two_d::fourier::chebyshev::transform <double> plan (horizontal_grid, vertical_grid, &data [0]);
	two_d::fourier::chebyshev::invert <double> iplan (horizontal_grid, vertical_grid, &data [0]);

	int flags;

	plan.execute (flags);	
	iplan.execute (flags);

	out.to_file ();

	
	// try {
	// 	element.run ();
	// } catch (...) {
	// 	FATAL ("Fatal error occurred. Check log.");
	// 	return 1;
	// }
	
	INFO ("Main complete.");
	
	return 0;
}
