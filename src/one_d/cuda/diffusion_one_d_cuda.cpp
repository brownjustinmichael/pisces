/*!**********************************************************************
 * \file diffusion_one_d_cuda.cpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-26.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "diffusion_one_d_cuda.hpp"

namespace cuda
{
	namespace one_d
	{
		template <class datatype>
		explicit_diffusion <datatype>::explicit_diffusion (datatype i_coeff, int i_n, bases::grid <datatype>* i_grid, datatype* i_data_in, datatype* i_data_out) :
		bases::explicit_plan <datatype> (i_n, i_data_in, i_data_out),
		coeff (i_coeff) {
		
			deriv_matrix.resize (n * n);
			deriv_matrix.copy_to_device (n * n, i_grid->get_data (2));
		
			TRACE ("Initialized.");
		}

		template <class datatype>
		void explicit_diffusion <datatype>::execute () {
			TRACE ("Operating...");
		
			// Set up and evaluate the explicit part of the diffusion equation
			utils::matrix_vector_multiply (n, n, coeff, deriv_matrix.pointer (), data_in, 1.0, data_out, n);

			TRACE ("Operation complete.");
		}
	
		template class explicit_diffusion <double>;
		template class explicit_diffusion <float>;
	} /* one_d */
} /* cuda */
