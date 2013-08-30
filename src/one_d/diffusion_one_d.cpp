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
	explicit_diffusion <datatype>::explicit_diffusion (datatype i_coeff, int i_n, bases::collocation_grid <datatype>* i_grid, datatype* i_data_in, datatype* i_data_out) :
	bases::explicit_plan <datatype> (i_n, i_data_in, i_data_out),
	coeff (i_coeff), 
	grid (i_grid) {
		TRACE ("Initialized.");
	}

	template <class datatype>
	void explicit_diffusion <datatype>::execute () {	
		TRACE ("Operating...");
		
		// if (*flags_ptr & linked_0) {
		// 	data_out [0] += coeff * utils::dot (n, grid->get_data (2), data_in, n);
		// }
		
		// Set up and evaluate the explicit part of the diffusion equation
		utils::matrix_vector_multiply (n, n, coeff, grid->get_data (2), data_in, 1.0, data_out, n);

		// if (*flags_ptr & linked_n) {
		// 	data_out [n - 1] += coeff * utils::dot (n, grid->get_data (2) + n - 1, data_in, n);
		// }

		TRACE ("Operation complete.");
	}
	
	template <class datatype>
	implicit_diffusion <datatype>::implicit_diffusion (datatype i_coeff, int i_n, bases::collocation_grid <datatype>* i_grid, datatype *i_matrix) : 
	bases::implicit_plan <datatype> (i_n, i_grid, i_matrix), 
	coeff (i_coeff) {}

	template <class datatype>
	void implicit_diffusion <datatype>::execute () {		
		TRACE ("Operating...");
	
		// This is the main loop for setting up the diffusion equation in Chebyshev space
		for (int i = 0; i < n; ++i) {
		   	utils::add_scaled (n, coeff, grid->get_data (2) + i, matrix + i, n, n);
		}
		TRACE ("Operation complete.");
	}
	
	template class explicit_diffusion <double>;
	template class explicit_diffusion <float>;

	template class implicit_diffusion <double>;
	template class implicit_diffusion <float>;
} /* one_d */

