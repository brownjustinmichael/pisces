/*!**********************************************************************
 * \file diffusion_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DIFFUSION_TWO_D_HPP_YHECX9VS
#define DIFFUSION_TWO_D_HPP_YHECX9VS

#include "../implicit_plan.hpp"

/*
	TODO Should non-axis diffusion be in the explicit or implicit rhs?
*/

namespace plans
{
	template <class datatype>
	class vertical_diffusion : public implicit_plan <datatype>
	{
	public:
		vertical_diffusion (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, datatype i_coeff, datatype i_alpha, datatype *i_matrix_n, datatype *i_matrix_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) : 
		implicit_plan <datatype> (i_grid_n, i_grid_m, i_matrix_n, i_matrix_m, i_data_in, i_data_out, i_element_flags, i_component_flags),
		coeff (i_coeff),
		alpha (i_alpha) {
			if (matrix_m) {
				for (int j = 0; j < m; ++j) {
					linalg::add_scaled (m, -coeff * alpha, grid_m.get_data (2) + j, matrix_m + j, m, m);
				}
			} else {
				WARN ("No matrix");
			}
		}
	
		virtual ~vertical_diffusion () {}

		void execute () {	
			TRACE ("Operating..." << element_flags);
			if (*component_flags & z_solve) {
				if (1.0 - alpha != 0.0) {
					linalg::matrix_matrix_multiply (m, ldn, m, coeff * (1.0 - alpha), grid_m.get_data (2), data_in, 1.0, data_out, m);
				}
			} else {
				linalg::matrix_matrix_multiply (m, ldn, m, coeff, grid_m.get_data (2), data_in, 1.0, data_out, m);
			}
			TRACE ("Operation complete.");
		}

		class factory : public implicit_plan <datatype>::factory
		{
		private:
			datatype coeff, alpha;
		public:
			factory (datatype i_coeff, datatype i_alpha) : coeff (i_coeff), alpha (i_alpha) {}
		
			virtual ~factory () {}
		
			virtual std::shared_ptr <plans::plan <datatype> > instance (plans::grid <datatype> **grids, datatype **matrices, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
				return std::shared_ptr <plans::plan <datatype> > (new vertical_diffusion <datatype> (*grids [0], *grids [1], coeff, alpha, matrices [0], matrices [1], i_data_in, i_data_out, i_element_flags, i_component_flags));
			}
		};

		using implicit_plan <datatype>::element_flags;
		using implicit_plan <datatype>::component_flags;

	private:
		datatype coeff;
		datatype alpha;
		using implicit_plan <datatype>::n;
		using implicit_plan <datatype>::ldn;
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
		finite_vertical_diffusion (plans::equation <datatype> &i_solver, datatype i_coeff, datatype i_alpha) :
		implicit_plan <datatype> (i_solver),
		coeff (i_coeff),
		alpha (i_alpha) {
			TRACE ("Instantiating...");
			for (int j = 0; j < m; ++j) {
				linalg::add_scaled (m, -coeff * alpha, grid_m.get_data (2) + j, matrix_m + j, m, m);
			}
		}
	
		virtual ~finite_vertical_diffusion () {}

		void execute () {	
			TRACE ("Operating...");
		
			if (*component_flags & z_solve) {
				if (1.0 - alpha != 0.0) {
					for (int i = 1; i < m - 1; ++i) {
						datatype pos_m1 = grid_m [i - 1];
						datatype pos_0 = grid_m [i];
						datatype pos_1 = grid_m [i + 1];
						for (int j = 0; j < ldn; ++j) {
							data_out [j * m + i] += coeff * (1.0 - alpha) * 2.0 * ((data_in [j * m + i + 1] - data_in [j * m + i]) / (pos_1 - pos_0) - (data_in [j * m + i] - data_in [j * m + i - 1]) / (pos_0 - pos_m1)) / (pos_1 - pos_m1);
						}
					}
				}
			} else {
				for (int i = 1; i < m - 1; ++i) {
					datatype pos_m1 = grid_m [i - 1];
					datatype pos_0 = grid_m [i];
					datatype pos_1 = grid_m [i + 1];
					for (int j = 0; j < ldn; ++j) {
						data_out [j * m + i] += coeff * 2.0 * ((data_in [j * m + i + 1] - data_in [j * m + i]) / (pos_1 - pos_0) - (data_in [j * m + i] - data_in [j * m + i - 1]) / (pos_0 - pos_m1)) / (pos_1 - pos_m1);
					}
				}
			}
		
			TRACE ("Operation complete.");
		}

		using implicit_plan <datatype>::element_flags;
		using implicit_plan <datatype>::component_flags;

	private:
		datatype coeff;
		datatype alpha;
		using implicit_plan <datatype>::n;
		using implicit_plan <datatype>::ldn;
		using implicit_plan <datatype>::m;
		using implicit_plan <datatype>::data_in;
		using implicit_plan <datatype>::data_out;
		using implicit_plan <datatype>::matrix_n;
		using implicit_plan <datatype>::matrix_m;
		using implicit_plan <datatype>::grid_n;
		using implicit_plan <datatype>::grid_m;
	};
} /* plans */

#endif /* end of include guard: DIFFUSION_TWO_D_HPP_YHECX9VS */
