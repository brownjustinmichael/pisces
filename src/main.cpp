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
#include "utils/profile.hpp"
#include "utils/rezone.hpp"
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
	int id = 0, n_elements = 1;

	// Initialize messenger
	bases::messenger process_messenger (&argc, &argv);

	try {
		id = process_messenger.get_id ();
		n_elements = process_messenger.get_np ();

		log_config::configure (&argc, &argv, id, "process_%d.log");
		std::string config_filename;
		
		if (argc <= 1) {
			config_filename = "../input/config.yaml";
		} else {
			config_filename = argv [1];
		}
		io::parameters config (config_filename);

		if (!config ["time.steps"].IsDefined ()) config ["time.steps"] = 1;
		if (!config ["grid.x.points"].IsDefined ()) config ["grid.x.points"] = 64;
		if (!config ["grid.z.points"].IsDefined ()) config ["grid.z.points"] = 64;
			
		omp_set_num_threads (config.get <int> ("parallel.maxthreads"));
		
		int m = config.get <int> ("grid.z.points") / n_elements + 1;
		
		std::vector <double> positions (n_elements + 1);
		for (int i = 0; i < n_elements + 1; ++i) {
			positions [i] = -config.get <double> ("grid.z.width") / 2.0 + config.get <double> ("grid.z.width") / n_elements * i;
		}
		
		int name = id;
		
		int n = config.get <int> ("grid.x.points");
		
		int n_steps = config.get <int> ("time.steps");
		
		bases::axis horizontal_axis (n, -config.get <double> ("grid.x.width") / 2.0, config.get <double> ("grid.x.width") / 2.0);
		bases::axis vertical_axis (m, positions [id], positions [id + 1], id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);
		
		// one_d::cosine::advection_diffusion_element <double> element (&vertical_axis, name, config, &process_messenger, 0x00);
		std::shared_ptr <bases::element <double>> element (new two_d::fourier::cosine::boussinesq_element <double> (horizontal_axis, vertical_axis, name, config, &process_messenger, 0x00));
		
		TRACE ("Element constructed.");
		
		if (config ["input.file"].IsDefined ()) {
			std::string file_format = "../input/" + config.get <std::string> ("input.file");
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);
			io::input input_stream (new io::two_d::netcdf (n, m), buffer);
		
			element->setup (&input_stream);
		}
		
		// Set up output
		std::shared_ptr <io::output> normal_stream;
		if (config ["output.file"].IsDefined ()) {
			std::string file_format = "../output/" + config.get <std::string> ("output.file");
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);
		
			normal_stream.reset (new io::incremental (new io::two_d::netcdf (n, m), buffer, config.get <int> ("output.every")));
		}
		
		std::shared_ptr <io::output> transform_stream;
		if (config ["output.transform_file"].IsDefined ()) {
			std::string file_format = "../output/" + config.get <std::string> ("output.transform_file");
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);
		
			transform_stream.reset (new io::incremental (new io::two_d::netcdf (n, m), buffer, config.get <int> ("output.every")));
		}
		
		/*
			TODO Because the output files are not within the element anymore, they become useless once the element is destroyed. Ponder this.
		*/
			
		clock_t cbegin, cend;
		std::chrono::time_point <std::chrono::system_clock> begin, end;
		
		cbegin = clock ();
		begin = std::chrono::system_clock::now ();

		while (n_steps > 0) {
			std::shared_ptr <io::virtual_dump> new_dump = element->rezone_minimize_ts (&positions [0], config.get <double> ("grid.rezone.min_size"), config.get <double> ("grid.rezone.max_size"), config.get <int> ("grid.rezone.n_tries"), config.get <int> ("grid.rezone.iters_fixed_t"), config.get <double> ("grid.rezone.step_size"), config.get <double> ("grid.rezone.k"), config.get <double> ("grid.rezone.t_initial"), config.get <double> ("grid.rezone.mu_t"), config.get <double> ("grid.rezone.t_min"));
			
			bases::axis vertical_axis (m, positions [id], positions [id + 1], id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);
			element.reset (new two_d::fourier::cosine::boussinesq_element <double> (horizontal_axis, vertical_axis, name, config, &process_messenger, 0x00));
			
			/*
				TODO It would be nice to combine the above construction of element with this one
			*/
			
			io::input *virtual_input (new io::input (new io::two_d::virtual_format (&*new_dump, n, m)));
			element->setup (&*virtual_input);
			if (normal_stream) {
				element->setup_output (normal_stream, normal_output);
			}
			if (transform_stream) {
				element->setup_output (transform_stream, transform_output);
			}
			
			// element.transform (forward_vertical | no_read);
			// for (int i = 0; i < config.get <int> ("time.steps"); ++i) {
			// 	element.transform (inverse_vertical | no_read | no_write);
			// 	element.transform (forward_vertical | no_read | no_write);
			// }
			// element.transform (inverse_vertical | no_write);

			element->run (n_steps, config.get <int> ("grid.rezone.check_every"));
		}
	
		cend = clock ();
		end = std::chrono::system_clock::now ();
			
		std::chrono::duration <double> eb = end - begin;
			
		INFO ("Main complete. CPU Time: " << ((double) (cend - cbegin))/CLOCKS_PER_SEC << " Wall Time: " << (double) eb.count () << " Efficiency: " << (((double) (cend - cbegin))/CLOCKS_PER_SEC / (double) eb.count () / omp_get_max_threads () * 100.) << "%");
	} catch (std::exception &except) {
		FATAL ("Fatal error occurred. Check log.");
		FATAL (except.what ());
		return 1;
		
		/*
			TODO Last check all should be somewhere not defined by the user
		*/
	} catch (int &except) {
		FATAL ("Fatal error occurred. Check log.");
		FATAL (except);
		return 1;
		
		/*
			TODO Last check all should be somewhere not defined by the user
		*/
	} catch (...) {
		FATAL ("Last ditch...");
		return 1;
	}
	
	return 0;
}
