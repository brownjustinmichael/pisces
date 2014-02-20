/*!**********************************************************************
 * \file main_cuda.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-08-14.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "config.hpp"
#include "one_d/cuda/element_one_d_cuda.hpp"
#include <vector>
#include <stdio.h>
#include "utils/solver_utils.hpp"

int main (int argc, char *argv[])
{	
	int id = 0, n_elements = 1;
	
	// Initialize messenger
	bases::messenger process_messenger (&argc, &argv, 2);
	
	id = process_messenger.get_id ();
	n_elements = process_messenger.get_np ();
	
	log_config::update_name (id);
	
	// The program runs through the execution flags.
	while ((argc > 1) && (argv [1] [0] == '-')) {
		switch (argv [1] [1]) {
			// Debug switch
			case 'D':
				log_config::update_severity (atoi (&(argv [1] [2])));
				break;
		}
		--argc;
		++argv;
	}
	
	TRACE ("Command line arguments read, beginning setup.");
	
		
	io::parameters params ("../input/parameters.txt");
	
	
	int n = params ["gridpoints"].asInt / n_elements;
	double position_0 = -1.0 + 2.0 / n_elements * id;
	double position_n = -1.0 + 2.0 / n_elements * (id + 1);
	int excess_0;
	int excess_n;
	if (id == 0) {
		excess_0 = 0;
	} else {
		excess_0 = 1;
	}
	if (id == n_elements - 1) {
		excess_n = 0;
	} else {
		excess_n = 1;
	}
	int name = id;
	
	if (id != 0) {
		TRACE ("Adding boundary to " << name << " at 0 at processor " << id - 1);
		process_messenger.add_boundary (one_d::edge_0, id - 1);
	}
	if (id != n_elements - 1) {
		TRACE ("Adding boundary to " << name << " at n - 1 at processor " << id + 1);
		process_messenger.add_boundary (one_d::edge_n, id + 1);
	}
	
	cuda::one_d::chebyshev::fft_element <double> element (n, excess_0, position_0, excess_n, position_n, name, params, &process_messenger, 0x00);
	
	element.setup ();
	
	try {
		element.run ();
	} catch (...) {
		FATAL ("Fatal error occurred. Check log.");
		return 0;
	}
	
	INFO ("Main complete.");
	
	return 0;
}