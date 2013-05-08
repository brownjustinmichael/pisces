/*!***********************************************************************
 * \file one_d/diffusion.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <iostream>
#include <cmath>
#include <cassert>
#include <memory>
#include "../bases/element.hpp"
#include "../bases/boundary.hpp"
#include "../config.hpp"
#include "diffusion_one_d.hpp"
#include "../utils/chebyshev.hpp"
#include "../utils/utils.hpp"

namespace one_d
{
	namespace chebyshev
	{
		explicit_diffusion::explicit_diffusion (double i_coeff, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, int i_name_in, int i_position, int i_name_out) : bases::explicit_plan (i_n, i_name_in, i_name_out) {
			MTRACE ("Initializing...");
			coeff = i_coeff;
			grid = i_grid;
			position = i_position;
			MTRACE ("Initialized.");
		}

		void explicit_diffusion::execute () {
			
			bases::explicit_plan::execute ();
		
			TRACE (logger, "Operating...");
			
			if (!(*flags_ptr & fixed_0)) {
				data_out [0] += coeff * utils::dot (n, grid->get_data (2), data_in, n);
			}
			
			// Set up and evaluate the explicit part of the diffusion equation
			utils::matrix_vector_multiply (n - 2, n, coeff, grid->get_data (2) + 1, data_in, 1.0, data_out + 1, n);

			if (!(*flags_ptr & fixed_n)) {
				data_out [n - 1] += coeff * utils::dot (n, grid->get_data (2) + n - 1, data_in, n);
			}

			TRACE (logger, "Operation complete.");
		}
		
		void explicit_diffusion::boundary (int index, bases::element* ext_element_ptr, int ext_index) {
			if (element_ptr && ext_element_ptr) {
				double x01 = ((&((*element_ptr) [position])) [index - 1] - (&((*element_ptr) [position])) [index]);
				double x02 = ((&((*element_ptr) [position])) [index - 1] - (&((*ext_element_ptr) [position])) [ext_index + 1]);
				double x12 = ((&((*element_ptr) [position])) [index] - (&((*ext_element_ptr) [position])) [ext_index + 1]);
				
				(&((*element_ptr) [name_out])) [index] += 2.0 * coeff * 2.0 * ((&((*element_ptr) [name_in])) [index - 1] * x12 - (&((*element_ptr) [name_in])) [index] * x02 + (&((*ext_element_ptr) [name_in])) [ext_index + 1] * x01) / x01 / x02 / x12;
				(&((*ext_element_ptr) [name_out])) [ext_index] += (&((*element_ptr) [name_out])) [index];
			}
		}
		

		implicit_diffusion::implicit_diffusion (double i_coeff, int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_matrix) : bases::implicit_plan (i_n, i_grid, i_matrix) {
			coeff = i_coeff;
			n = i_n;
			grid = i_grid;
			matrix = i_matrix;
		}

		void implicit_diffusion::execute () {
			if (!(*flags_ptr & factorized)) {
				TRACE (logger, "Operating...");
				bases::implicit_plan::execute ();
				// if (!(*flags_ptr & fixed_0)) {
				   	// utils::add_scaled (n, coeff, grid->get_data (2), matrix, n, n);
				// }
			
				// This is the main loop for setting up the diffusion equation in Chebyshev space
				for (int i = 1; i < n - 1; ++i) {
				   	utils::add_scaled (n, coeff, grid->get_data (2) + i, matrix + i, n, n);
				}
			
				// if (!(*flags_ptr & fixed_n)) {
				   	// utils::add_scaled (n, coeff, grid->get_data (2) + n - 1, matrix + n - 1, n, n);
				// }
				TRACE (logger, "Operation complete.");
			}
		}
	} /* chebyshev */
} /* one_d */

