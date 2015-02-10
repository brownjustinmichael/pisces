/*!**********************************************************************
 * \file config.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2015-01-27.
 * Copyright 2015 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef CONFIG_HPP_6DA8C3A3
#define CONFIG_HPP_6DA8C3A3

#include <string>
#include <map>
#include <vector>
#include <omp.h>
#include <algorithm>

#include "io/parameters.hpp"
#include "logger/logger.hpp"

io::parameters config (int *argc, char ***argv, int id, std::string override_config_filename = "config.yaml") {
	logger::log_config::configure (argc, argv, id, "process_%d.log");
	std::string config_filename;
	
	std::vector <std::string> config_strs;
	std::map <std::string, std::string> config_args;
	
	for (int i = 0; i < *argc; ++i) {
		if (((*argv) [i] [0] == '-') && ((*argv) [i] [1] == 'V')) {
			config_args [(*argv) [i + 1]] = (*argv) [i + 2];
			config_strs.push_back ((*argv) [i + 1]);
			(*argc) -= 3;
			for (int j = i; j < (*argc); ++j) {
				(*argv) [j] = (*argv) [j + 3];
			}
			i -= 1;
		}
	}
	
	if ((*argc) <= 1) {
		config_filename = override_config_filename;
	} else {
		config_filename = (*argv) [1];
	}
	INFO ("Reading from config file " << config_filename);
	io::parameters config (config_filename);
	
	for (int i = 0; i < (int) config_strs.size (); ++i) {
		INFO ("Overwriting config variable " << config_strs [i] << " with " << config_args [config_strs [i]]);
		config [config_strs [i]] = config_args [config_strs [i]];
	}
	
	omp_set_num_threads (std::min (config.get <int> ("parallel.maxthreads"), omp_get_max_threads ()));
	
	return config;
}

#endif /* end of include guard: CONFIG_HPP_6DA8C3A3 */
