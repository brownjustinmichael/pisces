/*!***********************************************************************
 * \file config.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "config.hpp"

int severity = 2;
log_config log_config_instance;

log_config::log_config () {}

void log_config::configure (int* argc, char*** argv, int id) {
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
