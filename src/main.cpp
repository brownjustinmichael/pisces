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
	bases::messenger process_messenger (&argc, &argv, 2);

	try {
		id = process_messenger.get_id ();
		n_elements = process_messenger.get_np ();

		log_config::configure (&argc, &argv, id);
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
		double position_m0 = -config.get <double> ("grid.z.width") / 2.0 + config.get <double> ("grid.z.width") / n_elements * id;
		double position_mm = -config.get <double> ("grid.z.width") / 2.0 + config.get <double> ("grid.z.width") / n_elements * (id + 1);
		int name = id;
		
		int n = config.get <int> ("grid.x.points");
		
		int n_steps = config.get <int> ("time.steps");
		
		bases::axis horizontal_axis (n, -config.get <double> ("grid.x.width") / 2.0, config.get <double> ("grid.x.width") / 2.0);
		bases::axis vertical_axis (m, position_m0, position_mm, id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);
		
		// one_d::cosine::advection_diffusion_element <double> element (&vertical_axis, name, config, &process_messenger, 0x00);
		std::shared_ptr <bases::element <double>> element (new two_d::fourier::cosine::boussinesq_element <double> (&horizontal_axis, &vertical_axis, name, config, &process_messenger, 0x00));
		
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
			element->setup_output (normal_stream, normal_output);
		}
		
		std::shared_ptr <io::output> transform_stream;
		if (config ["output.transform_file"].IsDefined ()) {
			std::string file_format = "../output/" + config.get <std::string> ("output.transform_file");
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);

			transform_stream.reset (new io::incremental (new io::two_d::netcdf (n, m), buffer, config.get <int> ("output.every")));
			element->setup_output (transform_stream, transform_output);
		}
		
		/*
			TODO Because the output files are not within the element anymore, they become useless once the element is destroyed. Ponder this.
		*/
	
		clock_t cbegin, cend;
		std::chrono::time_point <std::chrono::system_clock> begin, end;

		cbegin = clock ();
		begin = std::chrono::system_clock::now ();

		while (n_steps > 0) {
			try {
				element->run (n_steps);
			} catch (exceptions::mesh_adapt &except) {
				INFO ("Resetting grid...");
				element->write_transform_data ();
				element->transform (inverse_horizontal | inverse_vertical);
				element->read_transform_data ();
				/*
					TODO Remove these when updating to the more recent version
				*/
				
				io::virtual_dump dump;

				std::shared_ptr <io::output> virtual_output (new io::output (new io::two_d::virtual_format (&dump, n, m)));
				element->setup_output (virtual_output);
		
				virtual_output->to_file ();
				
				std::vector <double> position_profile (m);
				std::vector <double> velocity_profile (m);
				
				utils::profile <double> (m, n, &(dump.index <double> ("z")), &position_profile [0]);
				utils::profile <double> (m, n, &(dump.index <double> ("w")), &velocity_profile [0]);
				
				for (int i = 0; i < m; ++i) {
					velocity_profile [i]*= 100;
					DEBUG (velocity_profile [i]);
				}
				DEBUG (element->calculate_min_timestep (m, &position_profile [0], &velocity_profile [0], profile_timestep));

				io::input *virtual_input (new io::input (new io::two_d::virtual_format (&dump, n, m)));
				
				element.reset (new two_d::fourier::cosine::boussinesq_element <double> (&horizontal_axis, &vertical_axis, name, config, &process_messenger, 0x00));
				element->setup (&*virtual_input);
				if (normal_stream) {
					element->setup_output (normal_stream, normal_output);
				}
				if (transform_stream) {
					element->setup_output (transform_stream, transform_output);
				}
			}
		}
	
		cend = clock ();
		end = std::chrono::system_clock::now ();
	
		std::chrono::duration <double> eb = end - begin;
			
		INFO ("Main complete. CPU Time: " << ((double) (cend - cbegin))/CLOCKS_PER_SEC << " Wall Time: " << (double) eb.count () << " Efficiency: " << (((double) (cend - cbegin))/CLOCKS_PER_SEC / (double) eb.count () / omp_get_max_threads () * 100.) << "%");
	} catch (std::exception &except) {
		FATAL ("Fatal error occurred. Check log.");
		FATAL (except.what ());
		process_messenger.kill_all ();
		return 1;
		
		/*
			TODO Last check all should be somewhere not defined by the user
		*/
	} catch (int &except) {
		FATAL ("Fatal error occurred. Check log.");
		FATAL (except);
		process_messenger.kill_all ();
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
