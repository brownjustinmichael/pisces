/*!***********************************************************************
 * \file main.cpp
 * PISCES
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "mpi/messenger.hpp"

#include <memory>
#include <map>
#include <omp.h>
#include <ctime>
#include <chrono>
#ifdef VTRACE
#include "vt_user.h"
#endif

#include "logger/logger.hpp"
#include "io/parameters.hpp"
#include "plans/plan.hpp"
#include "elements/pseudo_element.hpp"
#include "elements/boussinesq.hpp"
#include "elements/rezone.hpp"
#include "config.hpp"

/*!*******************************************************************
 * \mainpage
 *
 * \author Justin Brown
 *
 * \section Introduction
 *
 * The goal of this project is to set up a code designed to do 2D 
 * Boussinesq simulations using a spectral element scheme.
 * 
 * WARNING: This code uses the C++ shared_ptr object. If you are unfamiliar with this object, this can lead to some problems. Always generate a shared_ptr either with the make_shared C++ standard library function or by setting the shared pointer to the result of the "new" operator. This will avoid many future headaches on your part. Never set a shared_ptr to a pointer of a previously existing object as the computer will attempt to call the destructor to the object twice. (This often appears as a segmentation fault after the code has finished while the computer cleans up.)
 *
 * Possible further reaching goals: 3D, pseudo-incompressible
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
		
		io::parameters parameters = config (&argc, &argv, id);
		
		int m = parameters.get <int> ("grid.z.points") / n_elements + 1;
		m += (m - 1) % 2;

		std::vector <double> positions (n_elements + 1);
		for (int i = 0; i < n_elements + 1; ++i) {
			positions [i] = -parameters.get <double> ("grid.z.width") / 2.0 + parameters.get <double> ("grid.z.width") / n_elements * i;
		}

		int name = id;

		int n = parameters.get <int> ("grid.x.points");

		grids::axis horizontal_axis (n, -parameters.get <double> ("grid.x.width") / 2.0, parameters.get <double> ("grid.x.width") / 2.0);
		grids::axis vertical_axis (m, positions [id], positions [id + 1], id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);
		
		TRACE ("Building data");
		
		data::thermo_compositional_data <double> data (&horizontal_axis, &vertical_axis, id, n_elements, parameters);
		
		TRACE ("Constructing element");
		
		std::shared_ptr <pisces::element <double>> element (new pisces::pseudo_element <double> (horizontal_axis, vertical_axis, name, parameters, data, &process_messenger, 0x00));
		
		if (pisces::element <double>::version () < versions::version ("0.6.0.0")) {
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
		std::shared_ptr <io::input> virtual_input;
		while (n_steps < parameters.get <int> ("time.steps") && element->duration < parameters.get <double> ("time.stop")) {
			if (parameters.get <int> ("grid.rezone.check_every") > 0 && n_steps != 0 && n_elements > 1) {
				INFO ("Rezoning");
				formats::virtual_file *virt = element->rezone_minimize_ts (&positions [0], parameters.get <double> ("grid.rezone.min_size"), parameters.get <double> ("grid.rezone.max_size"), parameters.get <int> ("grid.rezone.n_tries"), parameters.get <int> ("grid.rezone.iters_fixed_t"), parameters.get <double> ("grid.rezone.step_size"), parameters.get <double> ("grid.rezone.k"), parameters.get <double> ("grid.rezone.t_initial"), parameters.get <double> ("grid.rezone.mu_t"), parameters.get <double> ("grid.rezone.t_min"));
				
				if (virt) {
					formats::virtual_files ["main/virtual_file"] = *virt;
					grids::axis vertical_axis (m, positions [id], positions [id + 1], id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);
				
					virtual_input.reset (new io::formatted_input <formats::virtual_format> (formats::data_grid::two_d (n, m), "main/virtual_file"));
					data.setup (virtual_input);
				
					element.reset (new pisces::pseudo_element <double> (horizontal_axis, vertical_axis, name, parameters, data, &process_messenger, 0x00));
				}
			}
			element->run (n_steps, parameters.get <int> ("time.steps"), parameters.get <int> ("grid.rezone.check_every"), parameters.get <double> ("time.stop"));
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
