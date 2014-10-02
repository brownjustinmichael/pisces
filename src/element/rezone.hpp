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
#include "plan/grid.hpp"
#include "linalg/utils.hpp"

namespace utils
{
	template <class datatype>
	void rezone (utils::messenger *inter_messenger, bases::grid <datatype> *input_grid, bases::grid <datatype> *output_grid, io::formats::virtual_file *input_virtual_file, io::formats::virtual_file *output_virtual_file) {
		if (output_virtual_file != input_virtual_file) {
			*output_virtual_file = *input_virtual_file;
		}
		
		for (typename std::map <std::string, void *>::iterator iter = input_virtual_file->begin (); iter != input_virtual_file->end (); iter++) {
			if (input_virtual_file->dims [iter->first] [1] != 1 && input_virtual_file->check_type <datatype> (iter->first)) {
				TRACE ("Rezoning " << iter->first << "...");
				int nn = input_grid->get_n ();
				std::vector <int> ns (inter_messenger->get_np (), nn);
				inter_messenger->allgather <int> (1, &ns [0], &ns [0]);
				int nsum = 0;
				int nhere = 0;
				for (int i = 0; i < inter_messenger->get_np (); ++i) {
					nsum += ns [i];
				}
				std::vector <datatype> position_buffer (nsum, 0.0);
				std::vector <datatype> value_buffer (nsum * input_virtual_file->dims [iter->first] [0], 0.0);
				std::vector <datatype> inter_buffer (nsum * input_virtual_file->dims [iter->first] [0], 0.0);
			
				inter_messenger->allgatherv <datatype> (nn, &((*input_grid) [0]), &ns [0], &position_buffer [0]);
				for (int i = 0; i < inter_messenger->get_np (); ++i) {
					ns [i] *= input_virtual_file->dims [iter->first] [0];
					if (i < inter_messenger->get_id ()) {
						nhere += ns [i];
					}
				}
		
				linalg::matrix_switch (nn, input_virtual_file->dims [iter->first] [0], &(input_virtual_file->index <datatype> (iter->first)), &inter_buffer [nhere]);
				inter_messenger->allgatherv <datatype> (nn * input_virtual_file->dims [iter->first] [0], &inter_buffer [0], &ns [0], &inter_buffer [0]);
				linalg::matrix_switch (input_virtual_file->dims [iter->first] [0], nsum, &inter_buffer [0], &value_buffer [0]);
				output_virtual_file->add_var <datatype> (iter->first, input_virtual_file->dims [iter->first] [0], input_grid->get_n ());

				linalg::interpolate <datatype> (output_grid->get_n (), output_virtual_file->dims [iter->first] [0], nsum, 1.0, 0.0, &position_buffer [0], &value_buffer [0], &((*output_grid) [0]), &(output_virtual_file->index <datatype> (iter->first)));
				for (int i = 0; i < output_virtual_file->dims [iter->first] [0]; ++i) {
					for (int j = 0; j < output_virtual_file->dims [iter->first] [1]; ++j) {
						if (std::isnan ((&(output_virtual_file->index <datatype> (iter->first))) [i * output_virtual_file->dims [iter->first] [1] + j])) {
							FATAL ("Nan after interpolate");
							throw 0;
						}
					}
				}
			}
		}
	}
	
	template <class datatype>
	datatype minimum_timestep (int n, int m, bases::element <datatype> *element, utils::messenger *messenger, datatype *positions) {
		std::shared_ptr <io::output> virtual_output (new io::formatted_output <io::formats::two_d::virtual_format> (io::data_grid::two_d (n, m), "rezone/virtual_file", io::replace_file));
		element->setup_output (virtual_output);
		
		virtual_output->to_file ();
		
		int id = messenger->get_id ();
		int np = messenger->get_np ();
		
		bases::axis vertical_axis (m, positions [id], positions [id + 1], id == 0 ? 0 : 1, id == np - 1 ? 0 : 1);
		std::shared_ptr <bases::grid <double>> vertical_grid = element->generate_grid (&vertical_axis);
		
		rezone <datatype> (messenger, &*(element->grids [1]), &*vertical_grid, &io::virtual_files ["rezone/virtual_file"], &io::virtual_files ["rezone/new_virtual_file"]);
		
		return element->calculate_min_timestep (&io::virtual_files ["rezone/new_virtual_file"]);
	}
	
	
} /* utils */


#endif /* end of include guard: REZONE_CPP_UEW78T7Q */
