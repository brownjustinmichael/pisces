/*!***********************************************************************
 * \file master_one_d.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-05-14.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "../bases/master.hpp"

namespace one_d
{
	template <class Telement, class Tboundary>
	class master : public bases::master
	{
	public:
		master (int i_id, std::string parameter_filename, int i_n_elements, int* n_grid, double* position_grid, std::string* name_grid) : bases::master (i_id, i_n_elements, parameter_filename) {
			MTRACE ("Initializing...")
			for (int i = 0; i < i_n_elements; ++i) {
				MTRACE ("Adding element " << name_grid [i]);
				elements [i].reset (new Telement (n_grid [i] + 1, position_grid [i], position_grid [i + 1], name_grid [i], this->get_params (), 0x00));
				if (i != 0) {
					MTRACE ("Linking element " << name_grid [i - 1] << " at n - 1 with element " << name_grid [i] << " at 0");
					add_boundary (i - 1, linked_n, 2 * (i - 1) + 1, 2 * (i - 1) + 2);
					add_boundary (i, linked_0, 2 * (i - 1) + 2, 2 * (i - 1) + 1);
				}
			}
			MTRACE ("Initialized.");
		}
		virtual ~master () {}

		void add_boundary (int index, int edge, int send_id, int recv_id, int process = -1) {
			if (process == -1) {
				process = id;
			}
			elements [index]->add_boundary (std::make_shared <Tboundary> (Tboundary (edge, send_id, recv_id, process)));
		}
	};
} /* one_d */