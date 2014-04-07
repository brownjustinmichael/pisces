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

namespace utils
{
	template <class datatype>
	void rezone (bases::messenger *inter_messenger, bases::grid <datatype> *input_grid, bases::grid <datatype> *output_grid, io::virtual_dump *input_dump, io::virtual_dump *output_dump) {
		for (typename std::map <std::string, std::vector <datatype>>::iterator iter = input_dump->begin <datatype> (); iter != input_dump->end <datatype> (); iter++) {
			int nn = input_grid->n;
			std::vector <int> ns (inter_messenger->get_np (), nn);
			inter_messenger->allgather <int> (inter_messenger->get_np (), &ns [0]);
			std::vector <datatype> position_buffer (nn * input_dump->dims [iter->first] [1]);
			std::vector <datatype> value_buffer (nn * input_dump->dims [iter->first] [1]);
			inter_messenger->allgatherv <datatype> (nn, &((*input_grid) [0]), &ns [0], &position_buffer [0]);
			inter_messenger->allgatherv <datatype> (nn, &(input_dump->index <datatype> (iter->first)), &ns [0], &value_buffer [0]);
			output_dump->add_var <datatype> (iter->first, input_grid->n, input_dump->dims [iter->first] [1]);
			utils::interpolate (output_grid->n, output_dump->dims [iter->first] [1], nn, 1.0, &position_buffer [0], &value_buffer [0], &((*output_grid) [0]), &(output_dump->index <datatype> (iter->first, 0)));
		}
	}
} /* utils */


#endif /* end of include guard: REZONE_CPP_UEW78T7Q */
