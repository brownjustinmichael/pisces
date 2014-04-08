/*!**********************************************************************
 * \file rezone.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-25.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef REZONE_CPP_UEW78T7Q
#define REZONE_CPP_UEW78T7Q

#include "io.hpp"
#include "interpolate.hpp"
#include "../bases/grid.hpp"
#include "utils.hpp"

namespace utils
{
	template <class datatype>
	void rezone (bases::messenger *inter_messenger, bases::grid <datatype> *input_grid, bases::grid <datatype> *output_grid, io::virtual_dump *input_dump, io::virtual_dump *output_dump) {
		if (output_dump != input_dump) {
			*output_dump = *input_dump;
		}
		std::stringstream debug;
		
		for (typename std::map <std::string, std::vector <datatype>>::iterator iter = input_dump->begin <datatype> (); iter != input_dump->end <datatype> (); iter++) {
			if (input_dump->dims [iter->first] [1] != 1) {
				TRACE ("Rezoning " << iter->first << "...");
				int nn = input_grid->n;
				std::vector <int> ns (inter_messenger->get_np (), nn);
				inter_messenger->allgather <int> (1, &ns [0]);
				int nsum = 0;
				int nhere = 0;
				for (int i = 0; i < inter_messenger->get_np (); ++i) {
					nsum += ns [0];
				}
				std::vector <datatype> position_buffer (nsum);
				std::vector <datatype> value_buffer (nsum * input_dump->dims [iter->first] [0]);
				std::vector <datatype> inter_buffer (nsum * input_dump->dims [iter->first] [0]);
			
				inter_messenger->allgatherv <datatype> (nn, &((*input_grid) [0]), &ns [0], &position_buffer [0]);
				for (int i = 0; i < inter_messenger->get_np (); ++i) {
					ns [i] *= input_dump->dims [iter->first] [0];
					if (i < inter_messenger->get_id ()) {
						nhere += ns [i];
					}
				}

				utils::matrix_switch (nn, input_dump->dims [iter->first] [0], &(input_dump->index <datatype> (iter->first)), &inter_buffer [nhere]);
				inter_messenger->allgatherv <datatype> (nn * input_dump->dims [iter->first] [0], &inter_buffer [0], &ns [0]);
				utils::matrix_switch (input_dump->dims [iter->first] [0], nsum, &inter_buffer [0], &value_buffer [0]);
				output_dump->add_var <datatype> (iter->first, input_dump->dims [iter->first] [0], input_grid->n);
				for (int i = 0; i < nsum; ++i) {
					debug << position_buffer [i] << " ";
					for (int j = 0; j < input_dump->dims [iter->first] [0]; ++j) {
						debug << value_buffer [j * nsum + i] << " ";
					}
					DEBUG ("IN: " << debug.str ());
					debug.str ("");
				}
				DEBUG ("size " << output_dump->dims [iter->first] [0]);
				utils::interpolate (output_grid->n, output_dump->dims [iter->first] [0], nsum, 1.0, 0.0, &position_buffer [0], &value_buffer [0], &((*output_grid) [0]), &(output_dump->index <datatype> (iter->first)));
				for (int j = 0; j < output_grid->n; ++j) {
					debug << (*output_grid) [j] << " ";
					for (int i = 0; i < input_dump->dims [iter->first] [0]; ++i) {
						debug << output_dump->index <datatype> (iter->first, i, j) << " ";
					}
					DEBUG ("OUT: " << debug.str ());
					debug.str ("");
				}
			}
		}
	}
} /* utils */


#endif /* end of include guard: REZONE_CPP_UEW78T7Q */
