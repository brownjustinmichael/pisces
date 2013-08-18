/*!**********************************************************************
 * \file main_cuda.cpp
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-14.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "config.hpp"
#include "one_d/element_one_d.hpp"

#include "mpi.h"

int main (int argc, char *argv[])
{
	int id;
	int n_elements;
	
	// Initialize messenger
	bases::messenger process_messenger (&argc, &argv);

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
		
	io::parameter_map inputParams;
	io::read_params_txt parameters ("../input/parameters.txt");
	inputParams = parameters.load_params();
	
	
	int n = inputParams ["gridpoints"].asInt / n_elements;
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
	
	one_d::chebyshev::cuda_element element (n, position_0, position_n, excess_0, excess_n, name, inputParams, &process_messenger, 0x00);
	
	if (id != 0) {
		TRACE ("Adding boundary to " << name << " at 0 at processor " << id - 1);
		element.add_boundary (one_d::edge_0, 1, 2, id - 1);
	}
	if (id != n_elements - 1) {
		TRACE ("Adding boundary to " << name << " at n - 1 at processor " << id + 1);
		element.add_boundary (one_d::edge_n, 2, 1, id + 1);
	}
	
	element.send_positions ();
	
	element.run ();
	
	INFO ("Main complete.");
	
	return 0;
}