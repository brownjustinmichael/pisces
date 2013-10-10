/*!***********************************************************************
 * \file diffusion_one_d.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <iostream>
#include <cmath>
#include <cassert>
#include "../bases/element.hpp"
#include "../bases/collocation.hpp"
#include "../config.hpp"
#include "diffusion_one_d.hpp"
#include "../utils/utils.hpp"

namespace one_d
{
	template <class datatype>
	diffusion <datatype>::diffusion (int i_n, datatype i_coeff, datatype i_alpha, bases::collocation_grid <datatype>* i_grid, datatype* i_data_in, datatype* i_matrix, datatype* i_data_out, int i_flags) :
	bases::implicit_plan <datatype> (i_n, i_grid, i_data_in, i_matrix, i_data_out),
	coeff (i_coeff), 
	alpha (i_alpha),
	flags (i_flags) {
		for (int i = 0; i < n; ++i) {
			utils::add_scaled (n, -coeff * alpha, grid->get_data (2) + i, matrix + i, n, n);
		}
		flags |= implicit_set;
		TRACE ("Initialized.");
	}

	template <class datatype>
	void diffusion <datatype>::execute () {	
		TRACE ("Operating...");
		
		// Set up and evaluate the explicit part of the diffusion equation
		utils::matrix_vector_multiply (n, n, coeff * (1.0 - alpha), grid->get_data (2), data_in, 1.0, data_out, n);

		TRACE ("Operation complete.");
	}
	
	template class diffusion <double>;
	template class diffusion <float>;
} /* one_d */

