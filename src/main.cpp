/*!***********************************************************************
 * \file main.cpp
 * PISCES
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <memory>
#include <omp.h>
#include <ctime>
#include <chrono>
#ifdef VTRACE
#include "vt_user.h"
#endif

#include "logger/logger.hpp"
#include "mpi/messenger.hpp"
#include "io/parameters.hpp"
#include "io/input.hpp"
#include "io/output.hpp"
#include "io/formats/ascii.hpp"
#include "io/formats/netcdf.hpp"
#include "plans/plan.hpp"
#include "elements/boussinesq.hpp"
#include "elements/rezone.hpp"

/*!*******************************************************************
 * \mainpage
 *
 * \author Justin Brown
 * \author Ryan Moll
 *
 * \section Introduction
 *
 * The goal of this project is to set up a code designed to do 2D 
 * Anelastic simulations using a spectral element scheme.
 * 
 * WARNING: This code uses the C++ shared_ptr object. If you are unfamiliar with this object, this can lead to some problems. Always generate a shared_ptr either with the make_shared C++ standard library function or by setting the shared pointer to the result of the "new" operator. This will avoid many future headaches on your part. Never set a shared_ptr to a pointer of a previously existing object as the computer will attempt to call the destructor to the object twice. (This often appears as a segmentation fault after the code has finished while the computer cleans up.)
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
	mpi::messenger process_messenger (&argc, &argv);

	try {
		id = process_messenger.get_id ();
		n_elements = process_messenger.get_np ();

		logger::log_config::configure (&argc, &argv, id, "process_%d.log");
		std::string config_filename;

		if (argc <= 1) {
			config_filename = "config.yaml";
		} else {
			config_filename = argv [1];
		}
		io::parameters config (config_filename);

		omp_set_num_threads (config.get <int> ("parallel.maxthreads"));

		int m = config.get <int> ("grid.z.points") / n_elements + 1;
		m += m % 2;

		std::vector <double> positions (n_elements + 1);
		for (int i = 0; i < n_elements + 1; ++i) {
			positions [i] = -config.get <double> ("grid.z.width") / 2.0 + config.get <double> ("grid.z.width") / n_elements * i;
		}

		int name = id;

		int n = config.get <int> ("grid.x.points");

		plans::axis horizontal_axis (n, -config.get <double> ("grid.x.width") / 2.0, config.get <double> ("grid.x.width") / 2.0);
		plans::axis vertical_axis (m, positions [id], positions [id + 1], id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);

		data::thermo_compositional_data <double> data (&horizontal_axis, &vertical_axis, id, n_elements, config);
		
		std::shared_ptr <pisces::element <double>> element (new pisces::boussinesq_element <double> (horizontal_axis, vertical_axis, name, config, data, &process_messenger, 0x00));
		
		if (pisces::element <double>::version < versions::version ("0.6.0.0")) {
			INFO ("element.version < 0.6.0.0");
		}
		else {
			INFO ("element.version not < 0.6.0.0");
		}

		TRACE ("Element constructed.");

		clock_t cbegin, cend;
		std::chrono::time_point <std::chrono::system_clock> begin, end;

		cbegin = clock ();
		begin = std::chrono::system_clock::now ();

		int n_steps = 0;
		while (n_steps < config.get <int> ("time.steps")) {
			if (config.get <int> ("grid.rezone.check_every") > 0) {
				io::virtual_files ["main/virtual_file"] = *(element->rezone_minimize_ts (&positions [0], config.get <double> ("grid.rezone.min_size"), config.get <double> ("grid.rezone.max_size"), config.get <int> ("grid.rezone.n_tries"), config.get <int> ("grid.rezone.iters_fixed_t"), config.get <double> ("grid.rezone.step_size"), config.get <double> ("grid.rezone.k"), config.get <double> ("grid.rezone.t_initial"), config.get <double> ("grid.rezone.mu_t"), config.get <double> ("grid.rezone.t_min")));

				plans::axis vertical_axis (m, positions [id], positions [id + 1], id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);

				DEBUG (io::virtual_files ["main/virtual_file"].index <double> ("T", 16, 16));
				DEBUG (data.duration);
				DEBUG (io::virtual_files ["main/virtual_file"].index <double> ("t"));

				io::input *virtual_input (new io::formatted_input <io::formats::two_d::virtual_format> (io::data_grid::two_d (n, m), "main/virtual_file"));
				data.setup (&*virtual_input);
				DEBUG (data.duration);

				std::stringstream debug;
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j < m; ++j) {
						debug << *(data (temp, i, j)) << " ";
					}
					DEBUG ("OUT " << debug.str ());
					debug.str ("");
				}

				DEBUG ("DONE THERE");

				element.reset (new pisces::boussinesq_element <double> (horizontal_axis, vertical_axis, name, config, data, &process_messenger, 0x00));

				for (int j = 0; j < m; ++j) {
					for (int i = 0; i < n; ++i) {
						debug << data (x_velocity) [i * m + j] << " ";
					}
					DEBUG ("DATA U OUT " << debug.str ());
					debug.str ("");
				}

				/*
					TODO It would be nice to combine the above construction of element with this one
				*/
			}
			element->run (n_steps, config.get <int> ("time.steps"), config.get <int> ("grid.rezone.check_every"));
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
