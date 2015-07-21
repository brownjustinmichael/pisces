/*!***********************************************************************
 * \file element.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "element.hpp"

#include <cassert>
#include "logger/logger.hpp"
#include "linalg/utils.hpp"
#include <omp.h>
#include <yaml-cpp/yaml.h>
#include "linalg/exceptions.hpp"
#include "data/data.hpp"

#include <ctime>
#include <chrono>

/*!**********************************************************************
 * \brief Time a function call
 ************************************************************************/
#define TIME(call,cputime,duration) cbegin=clock(); tbegin=std::chrono::system_clock::now(); call; cend=clock(); tend=std::chrono::system_clock::now(); cputime+=((double) (cend - cbegin))/CLOCKS_PER_SEC; duration+=tend-tbegin;

namespace pisces
{
	template <class datatype>
	void element <datatype>::run (int &n_steps, int max_steps, int check_every, datatype stop) {
		TRACE ("Running...");
		datatype t_timestep;

		// Define the variables to use in timing the execution of the run
		clock_t cbegin, cend;
		std::chrono::time_point <std::chrono::system_clock> tbegin, tend;

		double transform_time = 0.0, execution_time = 0.0, solve_time = 0.0, factorize_time = 0.0, output_time = 0.0, timestep_time = 0.0;
		std::chrono::duration <double> transform_duration = std::chrono::duration <double>::zero (), execution_duration = std::chrono::duration <double>::zero (), solve_duration = std::chrono::duration <double>::zero (), factorize_duration = std::chrono::duration <double>::zero (), output_duration = std::chrono::duration <double>::zero (), timestep_duration = std::chrono::duration <double>::zero ();

		t_timestep = calculate_min_timestep ();
		messenger_ptr->min (&t_timestep);

		// Set up openmp to run multiple plans simultaneously
		omp_set_nested (true);
		int threads = params.get <int> ("parallel.maxthreads");
		std::stringstream debug;
		// Iterate through the total number of timesteps
		while (n_steps < max_steps && check_every != 0 && duration < stop) {
			data.reset ();
			INFO ("Timestep: " << n_steps);

			data.output ();
			
			// Factorize the matrices
			TIME (
			factorize ();
			, factorize_time, factorize_duration);

			TRACE ("Executing plans in spectral space...");

			TIME (
			for (iterator iter = begin (); iter != end (); iter++) {
				equations [*iter]->execute_plans (plans::solvers::pre_plan);
			}
			, execution_time, execution_duration);

			TRACE ("Executing plans in half-space...");

			TRACE ("Executing plans in real space...");

			// for (typename data::data <datatype>::iterator iter = data.begin (); iter < data.end (); ++iter)
			// {
			// 	if (data [*iter].update ()) {
			// 		if (data.transformers [*iter]) data.transformers [*iter]->update_from (real_real);
			// 	}
			// }
			
			TIME (
			t_timestep = calculate_min_timestep ();
			messenger_ptr->min (&t_timestep);
			, timestep_time, timestep_duration);

			TIME (
			for (iterator iter = begin (); iter != end (); iter++) {
				// DEBUG ("Executing");
				equations [*iter]->execute_plans (plans::solvers::post_plan);
			}
			, execution_time, execution_duration);

			TRACE ("Updating...");

			// Calculate the pre solver plans

			TIME (
			for (iterator iter = begin (); iter != end (); iter++) {
				equations [*iter]->execute_plans (plans::solvers::mid_plan);
			}
			, execution_time, execution_duration);

			data.output ();

			DEBUG ("WELL WELL BEFORE SOLVE " << (data ["temperature"].ptr (real_real) [200]));
			DEBUG ("WELL WELL BEFORE SOLVE " << (data ["temperature"].ptr (spectral_spectral) [200]));
			DEBUG ("WELL WELL BEFORE SOLVE " << (data ["temperature"].ptr (real_spectral) [200]));

			TIME (
			solve ();
			, solve_time, solve_duration);

			// Check whether the timestep has changed. If it has, mark all equations to be refactorized.

			duration += timestep;
			data.output ();

			INFO ("TOTAL TIME: " << duration);
			if (t_timestep != timestep) {
				for (std::vector <std::string>::iterator iter = data.begin (); iter != data.end (); iter++) {
					data [*iter].component_flags &= ~plans::solvers::factorized;
				}
				INFO ("Updating timestep: " << t_timestep);
			}
			timestep = t_timestep;

			TRACE ("Update complete");

			++n_steps;
			--check_every;
		}


		INFO ("Profiling Factorize: CPU Time: " << factorize_time << " Wall Time: " << (double) factorize_duration.count () << " Efficiency: " << factorize_time / (double) factorize_duration.count () / omp_get_max_threads () * 100. << "%");
		INFO ("Profiling Transform: CPU Time: " << transform_time << " Wall Time: " << (double) transform_duration.count () << " Efficiency: " << transform_time / (double) transform_duration.count () / omp_get_max_threads () * 100. << "%");
		INFO ("Profiling Execute: CPU Time: " << execution_time << " Wall Time: " << (double) execution_duration.count () << " Efficiency: " << execution_time / (double) execution_duration.count () / omp_get_max_threads () * 100. << "%");
		INFO ("Profiling Timestep: CPU Time: " << timestep_time << " Wall Time: " << (double) timestep_duration.count () << " Efficiency: " << timestep_time / (double) timestep_duration.count () / omp_get_max_threads () * 100. << "%");
		INFO ("Profiling Solve: CPU Time: " << solve_time << " Wall Time: " << (double) solve_duration.count () << " Efficiency: " << solve_time / (double) solve_duration.count () / omp_get_max_threads () * 100. << "%");
		INFO ("Profiling Output: CPU Time: " << output_time << " Wall Time: " << (double) output_duration.count () << " Efficiency: " << output_time / (double) output_duration.count () / omp_get_max_threads () * 100. << "%");
		INFO ("Max Threads = " << threads << " of " << omp_get_max_threads ());
	}
	
	template class element <double>;
} /* pisces */