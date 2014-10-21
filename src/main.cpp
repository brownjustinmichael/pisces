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

		if (!config ["time.steps"].IsDefined ()) config ["time.steps"] = 1;
		if (!config ["grid.x.points"].IsDefined ()) config ["grid.x.points"] = 64;
		if (!config ["grid.z.points"].IsDefined ()) config ["grid.z.points"] = 64;

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

		data::implemented_data <double> data (&horizontal_axis, &vertical_axis);
		
		std::shared_ptr <pisces::element <double>> element (new pisces::boussinesq_element <double> (horizontal_axis, vertical_axis, name, config, data, &process_messenger, 0x00));

		TRACE ("Element constructed.");

		const io::data_grid i_grid = io::data_grid::two_d (n, m, 0, config.get <bool> ("input.full") ? n_elements * m : 0, 0, config.get <bool> ("input.full") ? id * m : 0);
		
		if (config ["input.file"].IsDefined ()) {
			std::string file_format = "input/" + config.get <std::string> ("input.file");
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);
			io::formatted_input <io::formats::two_d::netcdf> input_stream (i_grid, buffer);

			data.setup (&input_stream);
		}
		
		const io::data_grid o_grid = io::data_grid::two_d (n, m, 0, config.get <bool> ("output.full") ? n_elements * m : 0, 0, config.get <bool> ("output.full") ? id * m : 0);
		
		// Set up output
		std::shared_ptr <io::output> normal_stream;
		if (config ["output.file"].IsDefined ()) {
			std::string file_format = "output/" + config.get <std::string> ("output.file");
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);

			normal_stream.reset (new io::appender_output <io::formats::two_d::netcdf> (o_grid, buffer, config.get <int> ("output.every")));
			data.setup_output (normal_stream);
		}

		std::shared_ptr <io::output> transform_stream;
		if (config ["output.transform_file"].IsDefined ()) {
			std::string file_format = "output/" + config.get <std::string> ("output.transform_file");
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);

			transform_stream.reset (new io::appender_output <io::formats::two_d::netcdf> (o_grid, buffer, config.get <int> ("output.every")));
			data.setup_output (transform_stream, transformed_horizontal | transformed_vertical);
		}

		std::shared_ptr <io::output> stat_stream;
		if (config ["output.stat.file"].IsDefined ()) {
			std::string file_format = "output/" + config.get <std::string> ("output.stat.file");
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);

			stat_stream.reset (new io::appender_output <io::formats::ascii> (io::data_grid::two_d (n, m), buffer, config.get <int> ("output.stat.every")));
			data.setup_stat (stat_stream);
		}

		/*
			TODO Setting up the streams should be more convenient
		*/

		clock_t cbegin, cend;
		std::chrono::time_point <std::chrono::system_clock> begin, end;

		cbegin = clock ();
		begin = std::chrono::system_clock::now ();

		int n_steps = 0;
		while (n_steps < config.get <int> ("time.steps")) {
			if (config.get <int> ("grid.rezone.check_every") > 0) {
				io::virtual_files ["main/virtual_file"] = *(element->rezone_minimize_ts (&positions [0], config.get <double> ("grid.rezone.min_size"), config.get <double> ("grid.rezone.max_size"), config.get <int> ("grid.rezone.n_tries"), config.get <int> ("grid.rezone.iters_fixed_t"), config.get <double> ("grid.rezone.step_size"), config.get <double> ("grid.rezone.k"), config.get <double> ("grid.rezone.t_initial"), config.get <double> ("grid.rezone.mu_t"), config.get <double> ("grid.rezone.t_min")));

				plans::axis vertical_axis (m, positions [id], positions [id + 1], id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);
				element.reset (new pisces::boussinesq_element <double> (horizontal_axis, vertical_axis, name, config, data, &process_messenger, 0x00));

				/*
					TODO It would be nice to combine the above construction of element with this one
				*/

				io::input *virtual_input (new io::formatted_input <io::formats::two_d::virtual_format> (io::data_grid::two_d (n, m), "main/virtual_file"));
				element->setup (&*virtual_input);
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
