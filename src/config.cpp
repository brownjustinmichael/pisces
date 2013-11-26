/*!***********************************************************************
 * \file config.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "config.hpp"

int severity = 2;

log_config::log_config (int* argc, char*** argv, int id) {
	while ((*argc > 1) && (*argv [1] [0] == '-')) {
		switch (*argv [1] [1]) {
			case 'D':
				severity = atoi (&(*argv [1] [2]));
				break;
		}
		--*argc;
		++*argv;
	}
}
