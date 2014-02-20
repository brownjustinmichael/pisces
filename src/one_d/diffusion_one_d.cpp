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
#include "../bases/grid.hpp"
#include "../config.hpp"
#include "diffusion_one_d.hpp"
#include "../utils/utils.hpp"

namespace one_d
{
	template <class datatype>
	diffusion <datatype>::diffusion (bases::grid <datatype> &i_grid, datatype i_coeff, datatype i_alpha, datatype *i_matrix, datatype *i_data_in, datatype *i_data_out, int *i_element_flags, int *i_component_flags) :
	implicit_plan <datatype> (i_grid, i_matrix, i_data_in, i_data_out, i_element_flags, i_component_flags),
	coeff (i_coeff), 
	alpha (i_alpha) {
		for (int i = 0; i < n; ++i) {
			utils::add_scaled (ld, -coeff * alpha, grid.get_data (2) + i, matrix + i, ld, n);
		}
		TRACE ("Initialized.");
	}

	template <class datatype>
	diffusion <datatype>::diffusion (bases::solver <datatype> &i_solver, datatype i_coeff, datatype i_alpha) :
	implicit_plan <datatype> (i_solver),
	coeff (i_coeff), 
	alpha (i_alpha) {
		for (int i = 0; i < n; ++i) {
			utils::add_scaled (ld, -coeff * alpha, grid.get_data (2) + i, matrix + i, ld, n);
		}
		TRACE ("Initialized.");
	}

	template <class datatype>
	void diffusion <datatype>::execute () {	
		TRACE ("Operating...");
		
		// Set up and evaluate the explicit part of the diffusion equation
		utils::matrix_vector_multiply (n, ld, coeff * (1.0 - alpha), grid.get_data (2), data_in, 1.0, data_out, ld);

		TRACE ("Operation complete.");
	}
	
	template class diffusion <double>;
	template class diffusion <float>;
	
	template <class datatype>
	nonlinear_diffusion <datatype>::nonlinear_diffusion (bases::solver <datatype> &i_solver, datatype i_coeff) :
	explicit_plan <datatype> (i_solver),
	coeff (i_coeff) {
		TRACE ("Initialized.");
	}

	template <class datatype>
	void nonlinear_diffusion <datatype>::execute () {	
		TRACE ("Operating...");
		
		// Set up and evaluate the explicit part of the nonlinear diffusion equation
		for (int i = 1; i < n - 1; ++i) {
			data_out [i] += 2.0 * coeff * data_in [i] * ((data_in [i + 1] - data_in [i]) / (grid [i + 1] - grid [i]) - (data_in [i] - data_in [i - 1]) / (grid [i] - grid [i - 1])) / (grid [i + 1] - grid [i - 1]);
			data_out [i] += coeff * (data_in [i + 1] - data_in [i - 1]) / (grid [i + 1] - grid [i - 1]) * (data_in [i + 1] - data_in [i - 1]) / (grid [i + 1] - grid [i - 1]);
		}

		TRACE ("Operation complete.");
	}
	
	template class nonlinear_diffusion <double>;
	template class nonlinear_diffusion <float>;
} /* one_d */

