/*!***********************************************************************
 * \file config.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "config.hpp"
#include <iostream>
#include <string>

int severity = 2;
log_config log_config_instance;

log_config::log_config () {}

void log_config::configure (int* argc, char*** argv, int id, std::string log_file) {
	for (int i = 0; i < *argc; ++i) {
		if (((*argv) [i] [0] == '-') && ((*argv) [i] [1] == 'D')) {
			severity = atoi (&((*argv) [i] [2]));
			--*argc;
			for (int j = i; j < *argc; ++j) {
				(*argv) [j] = (*argv) [j + 1];
			}
		}
	}
}

void log_config::trace (std::stringstream &log_stream) {
	std::string test = log_stream.str ();
	printf ("TRACE - %s\n", test.c_str ());
}

void log_config::debug (std::stringstream &log_stream) {
	std::string test = log_stream.str ();
	printf ("DEBUG - %s\n", test.c_str ());
}

void log_config::info (std::stringstream &log_stream) {
	std::string test = log_stream.str ();
	printf ("INFO - %s\n", test.c_str ());
}

void log_config::warn (std::stringstream &log_stream) {
	std::string test = log_stream.str ();
	printf ("WARN - %s\n", test.c_str ());
}

void log_config::error (std::stringstream &log_stream) {
	std::string test = log_stream.str ();
	printf ("ERROR - %s\n", test.c_str ());
}

void log_config::fatal (std::stringstream &log_stream) {
	std::string test = log_stream.str ();
	printf ("FATAL - %s\n", test.c_str ());
}