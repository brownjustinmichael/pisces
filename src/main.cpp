/*!***********************************************************************
 * \file main.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "bases/messenger.hpp"
#include "bases/plan.hpp"
#include "one_d/element_one_d.hpp"
#include "two_d/boussinesq_two_d.hpp"
#include "config.hpp"
#include "two_d/transform_two_d.hpp"
#include <memory>
#include <omp.h>
#include <ctime>
#include <chrono>
#ifdef VTRACE
#include "vt_user.h"
#endif

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
	try {
		int id = 0, n_elements = 1;

		// Initialize messenger
		bases::messenger process_messenger (&argc, &argv);

		id = process_messenger.get_id ();
		n_elements = process_messenger.get_np ();

		log_config::configure (&argc, &argv, id, "process_%d.log");
	
		io::parameters config ("../input/config.yaml");
		if (!config ["time.steps"].IsDefined ()) config ["time.steps"] = 1;
		if (!config ["grid.x.points"].IsDefined ()) config ["grid.x.points"] = 64;
		if (!config ["grid.z.points"].IsDefined ()) config ["grid.z.points"] = 64;
	
		omp_set_num_threads (config.get <int> ("parallel.maxthreads"));

		int m = config.get <int> ("grid.z.points") / n_elements + 1;
		double position_m0 = -config.get <double> ("grid.z.width") / 2.0 + config.get <double> ("grid.z.width") / n_elements * id;
		double position_mm = -config.get <double> ("grid.z.width") / 2.0 + config.get <double> ("grid.z.width") / n_elements * (id + 1);
		double position_n0 = -config.get <double> ("grid.x.width") / 2.0;
		double position_nn = config.get <double> ("grid.x.width") / 2.0;
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
		
		
		int n = config.get <int> ("grid.x.points");
	
		bases::axis horizontal_axis (n, position_n0, position_nn);
		bases::axis vertical_axis (m, position_m0, position_mm, excess_0, excess_n);
	
		one_d::cosine::advection_diffusion_element <double> element (&vertical_axis, name, config, &process_messenger, 0x00);
		// two_d::fourier::cosine::boussinesq_element <double> element (&horizontal_axis, &vertical_axis, name, config, &process_messenger, 0x00);

		clock_t cbegin, cend;
		std::chrono::time_point <std::chrono::system_clock> begin, end;

		cbegin = clock ();
		begin = std::chrono::system_clock::now ();

		element.transform (forward_vertical | no_read);
		for (int i = 0; i < config.get <int> ("time.steps"); ++i) {
			element.transform (inverse_vertical | no_read | no_write);
			element.transform (forward_vertical | no_read | no_write);
		}
		element.transform (inverse_vertical | no_write);
		
		// element.run (config.get <int> ("time.steps"));
	
		cend = clock ();
		end = std::chrono::system_clock::now ();
	
		std::chrono::duration <double> eb = end - begin;
	
		INFO ("Main complete. CPU Time: " << ((double) (cend - cbegin))/CLOCKS_PER_SEC << " Wall Time: " << (double) eb.count () << " Efficiency: " << (((double) (cend - cbegin))/CLOCKS_PER_SEC / (double) eb.count () / omp_get_max_threads () * 100.) << "%");
	} catch (std::exception& except) {
		FATAL (except.what ());
		FATAL ("Fatal error occurred. Check log.");
		return 1;
	}
	
	return 0;
}
