/*!**********************************************************************
 * \file rezone.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-25.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "rezone.hpp"

namespace pisces
{	
	void rezone (mpi::messenger *inter_messenger, grids::grid *input_grid, grids::grid *output_grid, formats::virtual_file *input_virtual_file, formats::virtual_file *output_virtual_file, double *value_buffer, double *inter_buffer) {
		if (output_virtual_file != input_virtual_file) {
			*output_virtual_file = *input_virtual_file;
		}
		
		int nn = input_grid->get_n ();
		std::vector <int> ns (inter_messenger->get_np (), nn);
		std::vector <int> nns (inter_messenger->get_np (), nn);
		inter_messenger->allgather <int> (1, &ns [0], &ns [0]);
		int nsum = 0;
		int nhere = 0;
		for (int i = 0; i < inter_messenger->get_np (); ++i) {
			nns [i] = ns [i];
			nsum += ns [i];
		}
		std::vector <double> position_buffer (nsum, 0.0);
		
		// Gather the global positions
		/*
			TODO Move this outside of the loop
		*/
		inter_messenger->allgatherv <double> (nn, &((*input_grid) [0]), &ns [0], &position_buffer [0]);
		
		// Iterate through the data
		for (std::map <std::string, void *>::iterator iter = input_virtual_file->begin (); iter != input_virtual_file->end (); iter++) {
			if (input_virtual_file->dims [iter->first] [1] != 1 && input_virtual_file->check_type <double> (iter->first)) {
				TRACE ("Rezoning " << iter->first << "...");
				
				nhere = 0;
				for (int i = 0; i < inter_messenger->get_np (); ++i) {
					nns [i] = ns [i] * input_virtual_file->dims [iter->first] [0];
					if (i < inter_messenger->get_id ()) {
						nhere += nns [i];
					}
				}
				
				// Gather the entirety of the values for the data in every element; need to switch to row-major order first
				linalg::matrix_switch (nn, input_virtual_file->dims [iter->first] [0], &(input_virtual_file->index <double> (iter->first)), &inter_buffer [nhere]);
				inter_messenger->allgatherv <double> (nn * input_virtual_file->dims [iter->first] [0], inter_buffer, &nns [0], inter_buffer);
				linalg::matrix_switch (input_virtual_file->dims [iter->first] [0], nsum, inter_buffer, value_buffer);
				output_virtual_file->add_var <double> (iter->first, input_virtual_file->dims [iter->first] [0], input_grid->get_n ());
				
				// Interpolate the new values from the global positions and values
				linalg::interpolate <double> (output_grid->get_n (), output_virtual_file->dims [iter->first] [0], nsum, 1.0, 0.0, &position_buffer [0], value_buffer, &((*output_grid) [0]), &(output_virtual_file->index <double> (iter->first)));
			}
		}
	}
} /* pisces */
