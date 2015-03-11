/*!**********************************************************************
 * \file rezone.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-25.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef REZONE_CPP_UEW78T7Q
#define REZONE_CPP_UEW78T7Q

#include "io/output.hpp"
#include "io/formats/virtual.hpp"
#include "linalg/interpolate.hpp"
#include "plans/grids/grid.hpp"
#include "linalg/utils.hpp"

namespace pisces
{
	/*!**********************************************************************
	 * \brief Rezone the elements according to the output grid
	 * 
	 * \param inter_messenger The mpi messenger that communicates between the elements
	 * \param input_grid A pointer to the grid describing the extent of the input data
	 * \param output_grid A pointer to the grid describing the desired extent of the output data
	 * \param input_virtual_file A pointer to the virtual file containing the input data
	 * \param output_virtual_file A pointer to the virtual file constructed by rezone. If NULL, use the input_virtual_file
	 * 
	 * This method takes an input_virtual_file and input_grid and rezones them according to the extent of output_grid. It will do so by communicating with the other elements to collect the data necessary for the rezone.
	 ************************************************************************/
	template <class datatype>
	void rezone (mpi::messenger *inter_messenger, grids::grid <datatype> *input_grid, grids::grid <datatype> *output_grid, formats::virtual_file *input_virtual_file, formats::virtual_file *output_virtual_file) {
		if (output_virtual_file != input_virtual_file) {
			*output_virtual_file = *input_virtual_file;
		}
		
		int nn = input_grid->get_n ();
		std::vector <int> ns (inter_messenger->get_np (), nn);
		inter_messenger->allgather <int> (1, &ns [0], &ns [0]);
		int nsum = 0;
		int nhere = 0;
		for (int i = 0; i < inter_messenger->get_np (); ++i) {
			nsum += ns [i];
		}
		std::vector <datatype> position_buffer (nsum, 0.0);
		
		// Gather the global positions
		/*
			TODO Move this outside of the loop
		*/
		inter_messenger->allgatherv <datatype> (nn, &((*input_grid) [0]), &ns [0], &position_buffer [0]);
		
		// Iterate through the data
		for (typename std::map <std::string, void *>::iterator iter = input_virtual_file->begin (); iter != input_virtual_file->end (); iter++) {
			if (input_virtual_file->dims [iter->first] [1] != 1 && input_virtual_file->check_type <datatype> (iter->first)) {
				TRACE ("Rezoning " << iter->first << "...");
				
				for (int i = 0; i < inter_messenger->get_np (); ++i) {
					ns [i] *= input_virtual_file->dims [iter->first] [0];
					if (i < inter_messenger->get_id ()) {
						nhere += ns [i];
					}
				}
				
				std::vector <datatype> value_buffer (nsum * input_virtual_file->dims [iter->first] [0], 0.0);
				std::vector <datatype> inter_buffer (nsum * input_virtual_file->dims [iter->first] [0], 0.0);
				
				// Gather the entirety of the values for the data in every element; need to switch to row-major order first
				linalg::matrix_switch (nn, input_virtual_file->dims [iter->first] [0], &(input_virtual_file->index <datatype> (iter->first)), &inter_buffer [nhere]);
				inter_messenger->allgatherv <datatype> (nn * input_virtual_file->dims [iter->first] [0], &inter_buffer [0], &ns [0], &inter_buffer [0]);
				linalg::matrix_switch (input_virtual_file->dims [iter->first] [0], nsum, &inter_buffer [0], &value_buffer [0]);
				output_virtual_file->add_var <datatype> (iter->first, input_virtual_file->dims [iter->first] [0], input_grid->get_n ());
				
				// Interpolate the new values from the global positions and values
				linalg::interpolate <datatype> (output_grid->get_n (), output_virtual_file->dims [iter->first] [0], nsum, 1.0, 0.0, &position_buffer [0], &value_buffer [0], &((*output_grid) [0]), &(output_virtual_file->index <datatype> (iter->first)));
			}
		}
	}
} /* pisces */


#endif /* end of include guard: REZONE_CPP_UEW78T7Q */
