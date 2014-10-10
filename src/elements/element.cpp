/*!***********************************************************************
 * \file bases/element.cpp
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

#include <ctime>
#include <chrono>

#ifndef VTRACE
#define VT_MARKER(x,y)
#define VT_MARKER_DEF(x,y) 0
#else
#include <vt_user.h>
#endif

/*!**********************************************************************
 * \brief Time a function call
 ************************************************************************/
#define TIME(call,cputime,duration) cbegin=clock(); tbegin=std::chrono::system_clock::now(); call; cend=clock(); tend=std::chrono::system_clock::now(); cputime+=((double) (cend - cbegin))/CLOCKS_PER_SEC; duration+=tend-tbegin;

namespace pisces
{
	template <class datatype>
	void element <datatype>::run (int &n_steps, int max_steps, int check_every) {
		TRACE ("Running...");
		datatype t_timestep;

		// Define the variables to use in timing the execution of the run
		clock_t cbegin, cend;
		std::chrono::time_point <std::chrono::system_clock> tbegin, tend;

		double transform_time = 0.0, execution_time = 0.0, solve_time = 0.0, factorize_time = 0.0, output_time = 0.0, timestep_time = 0.0;
		std::chrono::duration <double> transform_duration = std::chrono::duration <double>::zero (), execution_duration = std::chrono::duration <double>::zero (), solve_duration = std::chrono::duration <double>::zero (), factorize_duration = std::chrono::duration <double>::zero (), output_duration = std::chrono::duration <double>::zero (), timestep_duration = std::chrono::duration <double>::zero ();

		// If the element is in Cartesian space, transform to modal space and copy from the transform buffer
		transform (forward_horizontal | forward_vertical | no_read);

		// Set up openmp to run multiple plans simultaneously
		omp_set_nested (true);
		int threads = params.get <int> ("parallel.maxthreads");

		// Iterate through the total number of timesteps
		while (n_steps < max_steps && check_every != 0) {
			INFO ("Timestep: " << n_steps);

			// Transform the vertical grid to Cartesian space in the background
			TIME (
			transform (inverse_vertical | no_write | no_read | read_before);
			, transform_time, transform_duration);
	
			// Factorize the matrices
			TIME (
			factorize ();
			, factorize_time, factorize_duration);

			TRACE ("Executing plans...");

			TIME (
			for (iterator iter = begin (); iter != end (); iter++) {
				solvers [*iter]->execute_plans (pre_plan);
			}
			, execution_time, execution_duration);

			TIME (
			transform (inverse_horizontal | no_write | no_read | read_before);
			, transform_time, transform_duration);

			TIME (
			for (iterator iter = begin (); iter != end (); iter++) {
				solvers [*iter]->execute_plans (mid_plan);
			}
			, execution_time, execution_duration);
			
			TIME (
			transform (forward_horizontal | no_write | no_read | read_before);
			, transform_time, transform_duration);

			TIME (
			data.output ();
			data.output_stat ();
			for (int i = 0; i < (int) normal_profiles.size (); ++i) {
				normal_profiles [i]->to_file ();
			}
			, output_time, output_duration);

			// #pragma omp parallel sections num_threads(2)
				// {
				// #pragma omp section
					// {
			// Calculate the minimum timestep among all elements
						// omp_set_num_threads (threads);
						TIME (
						t_timestep = calculate_min_timestep ();
						messenger_ptr->min (&t_timestep);
						, timestep_time, timestep_duration);
					// }
				// #pragma omp section
					// {
						// omp_set_num_threads (threads);
						TIME (
						for (iterator iter = begin (); iter != end (); iter++) {
							solvers [*iter]->execute_plans (post_plan);
						}
						, execution_time, execution_duration);
					// }
				// }
			TRACE ("Updating...");

			// Transform forward in the horizontal direction
			TIME (
			transform (do_not_transform | no_write);
			, transform_time, transform_duration);

			// Calculate the pre solver plans
			TIME (
			for (iterator iter = begin (); iter != end (); iter++) {
				solvers [*iter]->execute_plans (pre_solve_plan);
			}
			, execution_time, execution_duration);

			// Output in transform space
			TIME (
			data.output_transform ();
			, output_time, output_duration);

			TIME (
			solve ();
			, solve_time, solve_duration);


			// Check whether the timestep has changed. If it has, mark all solvers to be refactorized.

			duration += timestep;
			INFO ("TOTAL TIME: " << duration);
			if (t_timestep != timestep) {
				for (std::map <int, int>::iterator iter = element_flags.begin (); iter != element_flags.end (); iter++) {
					element_flags [iter->first] &= ~factorized;
				}
				INFO ("Updating timestep: " << t_timestep);
			}
			timestep = t_timestep;

			TRACE ("Update complete");

			++n_steps;
			--check_every;
		}
		transform (do_not_transform | no_write);
		
		INFO ("Profiling Factorize: CPU Time: " << factorize_time << " Wall Time: " << (double) factorize_duration.count () << " Efficiency: " << factorize_time / (double) factorize_duration.count () / omp_get_max_threads () * 100. << "%");
		INFO ("Profiling Transform: CPU Time: " << transform_time << " Wall Time: " << (double) transform_duration.count () << " Efficiency: " << transform_time / (double) transform_duration.count () / omp_get_max_threads () * 100. << "%");
		INFO ("Profiling Execute: CPU Time: " << execution_time << " Wall Time: " << (double) execution_duration.count () << " Efficiency: " << execution_time / (double) execution_duration.count () / omp_get_max_threads () * 100. << "%");
		INFO ("Profiling Timestep: CPU Time: " << timestep_time << " Wall Time: " << (double) timestep_duration.count () << " Efficiency: " << timestep_time / (double) timestep_duration.count () / omp_get_max_threads () * 100. << "%");
		INFO ("Profiling Solve: CPU Time: " << solve_time << " Wall Time: " << (double) solve_duration.count () << " Efficiency: " << solve_time / (double) solve_duration.count () / omp_get_max_threads () * 100. << "%");
		INFO ("Profiling Output: CPU Time: " << output_time << " Wall Time: " << (double) output_duration.count () << " Efficiency: " << output_time / (double) output_duration.count () / omp_get_max_threads () * 100. << "%");
		INFO ("Max Threads = " << threads << " of " << omp_get_max_threads ());
	}
	
	template <class datatype>
	void element <datatype>::transform (int i_flags) {
		TRACE ("Transforming...");
		data.transform (i_flags);
// 		int threads = params.get <int> ("parallel.transform.threads");
// #pragma omp parallel num_threads (threads)
// 		{
// #pragma omp for
// 			for (int i = 0; i < (int) transforms.size (); ++i) {
// 				DEBUG ("Transforming " << transforms [i]);
// 				transformers [transforms [i]]->transform (i_flags);
// 			}
// 		}
	}
	
	template class element <double>;
} /* pisces */