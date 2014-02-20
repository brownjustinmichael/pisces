/*!***********************************************************************
 * \file bases/element.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "element.hpp"

#include <cassert>
#include "../config.hpp"
#include "solver.hpp"
#include "../utils/utils.hpp"
#include <omp.h>
#include <yaml-cpp/yaml.h>

#include <ctime>
#include <chrono>

#ifndef VTRACE
#define VT_MARKER(x,y)
#define VT_MARKER_DEF(x,y) 0
#else
#include <vt_user.h>
#endif

#define TIME(call,cputime,duration) cbegin=clock(); tbegin=std::chrono::system_clock::now(); call; cend=clock(); tend=std::chrono::system_clock::now(); cputime+=((double) (cend - cbegin))/CLOCKS_PER_SEC; duration+=tend-tbegin;

namespace bases
{
	template <class datatype>
	void element <datatype>::run (int n_steps) {
		TRACE ("Running...");
		datatype t_timestep;
		clock_t cbegin, cend;
		std::chrono::time_point <std::chrono::system_clock> tbegin, tend;
		
		double transform_time = 0.0, execution_time = 0.0, solve_time = 0.0, factorize_time = 0.0, output_time = 0.0, timestep_time = 0.0;
		std::chrono::duration <double> transform_duration, execution_duration, solve_duration, factorize_duration, output_duration, timestep_duration;
	
		write_transform_data ();
		transform (forward_horizontal | forward_vertical);
		
		omp_set_nested (true);
		int threads = params.get <int> ("parallel.maxthreads");
		
		for (int j = 0; j <= n_steps; ++j) {
			TIME (
			read_transform_data ();
			transform (inverse_vertical);
			, transform_time, transform_duration);

			TIME (
			factorize ();
			, factorize_time, factorize_duration);
			
			INFO ("Timestep " << j << " of " << n_steps);
	
			TRACE ("Calculating...");
	
			explicit_reset ();
			
			TRACE ("Executing plans...");

			TIME (
			for (iterator iter = begin (); iter != end (); iter++) {
				solvers [*iter]->execute_pre_plans ();
			}
			, execution_time, execution_duration);
			
			TIME (
			read_transform_data ();
			transform (inverse_horizontal);
			, transform_time, transform_duration);

			TIME (
			for (iterator iter = begin (); iter != end (); iter++) {
				solvers [*iter]->execute_mid_plans ();
			}
			, execution_time, execution_duration);
	
			TIME (
			read_transform_data ();
			transform (forward_horizontal);
			, transform_time, transform_duration);

			if (normal_stream) {
			#pragma omp parallel sections num_threads(3)
				{
				#pragma omp section
					{
			// Calculate the minimum timestep among all elements
						omp_set_num_threads (threads);
						TIME (
						t_timestep = calculate_min_timestep ();
						messenger_ptr->min (&t_timestep);
						, timestep_time, timestep_duration);
					}
				#pragma omp section
					{
						TIME (
						if (normal_stream) {
							normal_stream->to_file ();
						}
						, output_time, output_duration);
					}
				#pragma omp section
					{
						omp_set_num_threads (threads);
						TIME (
						for (iterator iter = begin (); iter != end (); iter++) {
							solvers [*iter]->execute_post_plans ();
						}
						, execution_time, execution_duration);
						
					}
				}
			} else {
			#pragma omp parallel sections num_threads(2)
				{
				#pragma omp section
					{
			// Calculate the minimum timestep among all elements
						omp_set_num_threads (threads);
						TIME (
						t_timestep = calculate_min_timestep ();
						messenger_ptr->min (&t_timestep);
						, timestep_time, timestep_duration);
					}
				#pragma omp section
					{
						omp_set_num_threads (threads);
						TIME (
						for (iterator iter = begin (); iter != end (); iter++) {
							solvers [*iter]->execute_post_plans ();
						}
						, execution_time, execution_duration);
					}
				}
			}

			TRACE ("Updating...");

			// Transform forward in the horizontal direction
			TIME (
			read_transform_data ();
			, transform_time, transform_duration);
	
			// Calculate the pre solver plans
			TIME (
			for (iterator iter = begin (); iter != end (); iter++) {
				solvers [*iter]->execute_pre_solve_plans ();
			}
			, execution_time, execution_duration);
			
			// Output in transform space
			TIME (
			if (transform_stream) {
				TRACE ("Writing to file...");
				transform_stream->to_file ();
			}
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
		}
		
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
		typedef typename std::vector <int>::iterator iterator;
		omp_set_nested (true);
		
		int threads = params.get <int> ("parallel.transform.threads");
		
		#pragma omp parallel for num_threads (threads)
		for (int i = 0; i < (int) transforms.size (); ++i) {
			master_transforms [transforms [i]]->transform (i_flags);
		}
	}
	
	template <class datatype>
	void element <datatype>::write_transform_data () {
		TRACE ("Writing to GPU...");
		typedef typename std::vector <int>::iterator iterator;
		omp_set_nested (true);
		
		int threads = params.get <int> ("parallel.transform.threads");
		
		#pragma omp parallel for num_threads (threads)
		for (int i = 0; i < (int) transforms.size (); ++i) {
			master_transforms [transforms [i]]->write ();
		}
	}
	
	template <class datatype>
	void element <datatype>::read_transform_data () {
		TRACE ("Reading from GPU...");
		typedef typename std::vector <int>::iterator iterator;
		omp_set_nested (true);
		
		int threads = params.get <int> ("parallel.transform.threads");
		
		#pragma omp parallel for num_threads (threads)
		for (int i = 0; i < (int) transforms.size (); ++i) {
			master_transforms [transforms [i]]->read ();
		}
	}
	
	template class element <double>;
} /* bases */