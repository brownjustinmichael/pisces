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

/*
	TODO Should non-axis diffusion be in the explicit or implicit rhs?
*/

namespace two_d
{
	namespace fourier
	{
		namespace chebyshev
		{
			template <class datatype>
			class horizontal_diffusion : public implicit_plan <datatype>
			{
			public:
				horizontal_diffusion (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype i_coeff, datatype i_alpha, datatype *i_matrix_n, datatype *i_matrix_m, datatype* i_data_in, datatype* i_data_out = NULL, int i_flags = 0x0) :
				implicit_plan <datatype> (i_grid_n, i_grid_m, i_matrix_n, i_matrix_m, i_data_in, i_data_out),
				coeff (i_coeff),
				alpha (i_alpha),
				flags (i_flags) {
					pioL2 = 4.0 * (std::acos (-1.0) * std::acos (-1.0) / (i_grid_n.position (n - 1) - i_grid_n.position (0)) / (i_grid_n.position (n - 1) - i_grid_n.position (0)));
					for (int i = 0; i < 2 * (n / 2 + 1); ++i) {
						matrix_n [i] = coeff * alpha * pioL2 * (datatype) ((i / 2) * (i / 2));
					}
				}
				
				virtual ~horizontal_diffusion () {}
			
				void execute (int &element_flags) {	
					TRACE ("Operating...");
					std::stringstream debug;
					if (element_flags & x_solve) {
						if (1.0 - alpha != 0.0) {
							for (int j = 0; j < m; ++j) {
								for (int i = 0; i < 2 * (n / 2 + 1); ++i) {
									data_out [i * m + j] -= coeff * (1.0 - alpha) * pioL2 * (datatype) ((i / 2) * (i / 2)) * data_in [i * m + j];
								}
							}
						}
					} else {
						for (int j = 0; j < m; ++j) {
							for (int i = 0; i < 2 * (n / 2 + 1); ++i) {
								data_out [i * m + j] -= coeff * pioL2 * (datatype) ((i / 2) * (i / 2)) * data_in [i * m + j];
							}
						}
					}
					TRACE ("Operation complete.");
				}
			
			private:
				datatype coeff;
				datatype alpha;
				datatype pioL2;
				int flags;
				using implicit_plan <datatype>::n;
				using implicit_plan <datatype>::m;
				using implicit_plan <datatype>::data_in;
				using implicit_plan <datatype>::data_out;
				using implicit_plan <datatype>::matrix_n;
			};
			
			template <class datatype>
			class vertical_diffusion : public implicit_plan <datatype>
			{
			public:
				vertical_diffusion (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype i_coeff, datatype i_alpha, datatype *i_matrix_n, datatype *i_matrix_m, datatype* i_data_in, datatype* i_data_out = NULL, int i_flags = 0x0) :
				implicit_plan <datatype> (i_grid_n, i_grid_m, i_matrix_n, i_matrix_m, i_data_in, i_data_out),
				coeff (i_coeff),
				alpha (i_alpha),
				flags (i_flags) {
					for (int j = 0; j < m; ++j) {
						utils::add_scaled (m, -coeff * alpha, grid_m.get_data (2) + j, matrix_m + j, m, m);
					}
				}
				
				virtual ~vertical_diffusion () {}
			
				void execute (int &element_flags) {	
					TRACE ("Operating...");
					
					if (element_flags & z_solve) {
						if (1.0 - alpha != 0.0) {
							utils::matrix_matrix_multiply (m, 2 * (n / 2 + 1), m, coeff * (1.0 - alpha), grid_m.get_data (2), data_in, 1.0, data_out, m);
						}
					} else {
						utils::matrix_matrix_multiply (m, 2 * (n / 2 + 1), m, coeff, grid_m.get_data (2), data_in, 1.0, data_out, m);
					}
					
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
			
			template <class datatype>
			class finite_vertical_diffusion : public implicit_plan <datatype>
			{
			public:
				finite_vertical_diffusion (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype i_coeff, datatype i_alpha, datatype *i_matrix_n, datatype *i_matrix_m, datatype* i_data_in, datatype* i_data_out = NULL, int i_flags = 0x0) :
				implicit_plan <datatype> (i_grid_n, i_grid_m, i_matrix_n, i_matrix_m, i_data_in, i_data_out),
				coeff (i_coeff),
				alpha (i_alpha),
				flags (i_flags) {
					for (int j = 0; j < m; ++j) {
						utils::add_scaled (m, -coeff * alpha, grid_m.get_data (2) + j, matrix_m + j, m, m);
					}
				}
				
				virtual ~finite_vertical_diffusion () {}
			
				void execute (int &element_flags) {	
					TRACE ("Operating...");
					
					if (element_flags & z_solve) {
						if (1.0 - alpha != 0.0) {
							for (int i = 1; i < m - 1; ++i) {
								datatype pos_m1 = grid_m.position (i - 1);
								datatype pos_0 = grid_m.position (i);
								datatype pos_1 = grid_m.position (i + 1);
								for (int j = 0; j < 2 * (n / 2 + 1); ++j) {
									data_out [j * m + i] += coeff * (1.0 - alpha) * 2.0 * ((data_in [j * m + i + 1] - data_in [j * m + i]) / (pos_1 - pos_0) - (data_in [j * m + i] - data_in [j * m + i - 1]) / (pos_0 - pos_m1)) / (pos_1 - pos_m1);
								}
							}
						}
					} else {
						for (int i = 1; i < m - 1; ++i) {
							datatype pos_m1 = grid_m.position (i - 1);
							datatype pos_0 = grid_m.position (i);
							datatype pos_1 = grid_m.position (i + 1);
							for (int j = 0; j < 2 * (n / 2 + 1); ++j) {
								data_out [j * m + i] += coeff * 2.0 * ((data_in [j * m + i + 1] - data_in [j * m + i]) / (pos_1 - pos_0) - (data_in [j * m + i] - data_in [j * m + i - 1]) / (pos_0 - pos_m1)) / (pos_1 - pos_m1);
							}
						}
					}
					
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
