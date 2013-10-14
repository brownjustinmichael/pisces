/*!**********************************************************************
 * \file diffusion_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DIFFUSION_TWO_D_HPP_YHECX9VS
#define DIFFUSION_TWO_D_HPP_YHECX9VS

#include "../bases/plan.hpp"
#include "solver_two_d.hpp"

namespace two_d
{
	namespace fourier
	{
		namespace chebyshev
		{
			template <class datatype>
			class diffusion : public implicit_plan
			{
			public:
				diffusion (int i_n, int i_m, datatype i_coeff, datatype i_alpha, bases::grid <datatype>* i_grid_n, bases::grid <datatype>* i_grid_m, datatype* i_data_in, datatype* i_matrix_n, datatype* i_matrix_m, datatype* i_data_out = NULL, int i_flags = 0x0) :
				implicit_plan <datatype> (i_n, i_m, i_grid_n, i_grid_m, i_data_in, i_matrix_n, i_matrix_m, i_data_out),
				coeff (i_coeff),
				alpha (i_alpha),
				flags (i_flags) {
					datatype pioM = -coeff * alpha * (2.0 * std::acos (-1.0) / (n - 1)) ** 2;
					for (int i = 0; i < n; ++j) {
						matrix_n [j] = pioM * (datatype) ((i / 2) * (i / 2));
					}
					for (int j = 0; j < m; ++j) {
						utils::add_scaled (m, -coeff * alpha, grid_m->get_data (2) + i, matrix_m + i, m, m);
					}
				}
				
				virtual ~diffusion () {}
			
				execute () {	
					TRACE ("Operating...");
					
					if (flags & x_solve) {
						utils::matrix_vector_multiply (n, n, coeff * (1.0 - alpha), grid_n->get_data (2), data_in, 1.0, data_out, n);
					} else {
						utils::matrix_vector_multiply (n, n, coeff, grid_n->get_data (2), data_in, 1.0, data_out, n);
					}
					
					if (flags & y_solve) {
						utils::matrix_vector_multiply (m, m, coeff * (1.0 - alpha), grid_m->get_data (2), data_in, 1.0, data_out, m);
					} else {
						utils::matrix_vector_multiply (m, m, coeff, grid_m->get_data (2), data_in, 1.0, data_out, m);
					}
					
					/*
						TODO Make more efficient by doing calculation in spectral space, at least for horizontal component
					*/

					TRACE ("Operation complete.");
				}
			
			private:
				datatype coeff;
				datatype alpha;
				int flags;
				using implicit_plan <datatype>::n;
				using implicit_plan <datatype>::m;
				using implicit_plan <datatype>::data_in;
				using implicit_plan <datatype>::data_out;
				using implicit_plan <datatype>::matrix_n;
				using implicit_plan <datatype>::matrix_m;
				using implicit_plan <datatype>::grid_n;
				using implicit_plan <datatype>::grid_m;
			};
		} /* chebyshev */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: DIFFUSION_TWO_D_HPP_YHECX9VS */
