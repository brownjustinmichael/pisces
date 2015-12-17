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
#include "rezone.hpp"

#include <ctime>
#include <chrono>

/*!**********************************************************************
 * \brief Time a function call
 ************************************************************************/
#define TIME(call,cputime,duration) cbegin=clock(); tbegin=std::chrono::system_clock::now(); call; cend=clock(); tend=std::chrono::system_clock::now(); cputime+=((double) (cend - cbegin))/CLOCKS_PER_SEC; duration+=tend-tbegin;

namespace pisces
{
	void element::run (int &n_steps) {
		TRACE ("Running...");
		double t_timestep;

		// Define the variables to use in timing the execution of the run
		clock_t cbegin, cend;
		std::chrono::time_point <std::chrono::system_clock> tbegin, tend;

		double transform_time = 0.0, execution_time = 0.0, solve_time = 0.0, factorize_time = 0.0, output_time = 0.0, timestep_time = 0.0;
		std::chrono::duration <double> transform_duration = std::chrono::duration <double>::zero (), execution_duration = std::chrono::duration <double>::zero (), solve_duration = std::chrono::duration <double>::zero (), factorize_duration = std::chrono::duration <double>::zero (), output_duration = std::chrono::duration <double>::zero (), timestep_duration = std::chrono::duration <double>::zero ();

		data.add_stream_attribute("version", version());
		data.add_stream_attribute("element", class_name());
		data.add_stream_attribute("pid", std::to_string ((long long int) messenger_ptr->get_id()));
		data.add_stream_attribute("np", std::to_string ((long long int) messenger_ptr->get_np()));

		t_timestep = calculate_min_timestep ();
		messenger_ptr->min (&t_timestep);

		for (data::data::iterator iter = data.begin (); iter != data.end (); ++iter) {
			if (data.transformers [*iter]) data.transformers [*iter]->update ();
		}

		// Set up openmp to run multiple plans simultaneously
		omp_set_nested (true);
		int threads = params.get <int> ("parallel.maxthreads");
		std::stringstream debug;
		double stop = params ["time.stop"].as <double> ();
		int max_steps = params.get <int> ("time.steps");
		int check_every = params.get <int> ("grid.rezone.check_every");
		// Iterate through the total number of timesteps
		while (n_steps < max_steps && check_every != 0 && duration < stop) {
			data.n_steps = n_steps;
			data.reset ();
			INFO ("Timestep: " << n_steps);

			TIME (
			data.output ();
			, output_time, output_duration);
			
			TRACE ("Executing plans...");

			TIME (
			for (iterator iter = begin (); iter != end (); iter++) {
				equations [*iter]->execute_plans ();
			}
			, execution_time, execution_duration);

			TIME (
			t_timestep = calculate_min_timestep ();
			messenger_ptr->min (&t_timestep);
			, timestep_time, timestep_duration);

			// Factorize the matrices
			TIME (
			factorize ();
			, factorize_time, factorize_duration);

			INFO ("TOTAL TIME: " << duration);
			if (t_timestep + duration > stop) t_timestep = stop - duration;
			if (t_timestep != timestep) {
				for (std::vector <std::string>::iterator iter = data.begin (); iter != data.end (); iter++) {
					data [*iter].component_flags &= ~plans::solvers::factorized;
				}
				INFO ("Updating timestep: " << t_timestep);
			}
			timestep = t_timestep;

			TIME (
			solve ();
			, solve_time, solve_duration);

			// Check whether the timestep has changed. If it has, mark all equations to be refactorized.

			duration += timestep;

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

	formats::virtual_file *element::rezone_minimize_ts (double *positions, double min_size, double max_size, int n_tries, int iters_fixed_t, double step_size, double k, double t_initial, double mu_t, double t_min, double (*function) (element *, formats::virtual_file *)) {
		TRACE ("Rezoning...");

		rezone_virtual_file = make_virtual_file (profile_only | timestep_only);
		
		rezone_data data;
		data.element_ptr = this;
		data.id = messenger_ptr->get_id ();
		data.np = messenger_ptr->get_np ();
		data.merit_func = function;
		data.min_size = min_size;
		data.max_size = max_size;
		
		get_zoning_positions (data.positions);
		
		double original = rezone_data::func (&data);
		
		gsl_siman_params_t params = {n_tries, iters_fixed_t, step_size, k, t_initial, mu_t, t_min};

		const gsl_rng_type * T;
		gsl_rng * r;
		
		gsl_rng_env_setup();
		
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);

		gsl_siman_solve(r, &data, rezone_data::func, rezone_data::rezone_generate_step, rezone_data::rezone_step_size, NULL, NULL, NULL, NULL, sizeof (rezone_data), params);

		gsl_rng_free (r);
		
		double value = rezone_data::func (&data);
		
		DEBUG ("Old: " << original << " New: " << value);
		
		if (value < original * (original < 0. ? 1. / rezone_mult : rezone_mult)) {
			for (int i = 0; i < data.np + 1; ++i) {
				positions [i] = data.positions [i];
			}
			return make_rezoned_virtual_file (data.positions, make_virtual_file ());
		}
		
		return NULL;
	}
} /* pisces */