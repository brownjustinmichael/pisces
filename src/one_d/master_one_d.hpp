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
			double diffusion_coeff = this->get_params () ["diffusion_coeff"].asDouble;
			
			if (i_n_elements > 1) {
				buffers1.resize (i_n_elements - 1);
				buffers2.resize (i_n_elements - 1);
			}
			for (int i = 0; i < i_n_elements; ++i) {				
				elements [i].reset (new Telement (n_grid [i] + 1, position_grid [i], position_grid [i + 1], name_grid [i], this->get_params (), 0x00));
				if (i != 0) {
					buffers1 [i - 1].resize (100);
					buffers2 [i - 1].resize (100);
					elements [i - 1]->add_boundary (std::make_shared <Tboundary> (Tboundary (linked_0, diffusion_coeff, position, velocity, rhs, &(buffers1 [i - 1] [0]), &(buffers2 [i - 1] [0]))));
					elements [i]->add_boundary (std::make_shared <Tboundary> (Tboundary (linked_n, diffusion_coeff, position, velocity, rhs, &(buffers2 [i - 1] [0]), &(buffers1 [i - 1] [0]))));
				}
			}
		}
		virtual ~master () {}
	
	private:
		std::vector <std::vector <double>> buffers1, buffers2;
	};
} /* one_d */