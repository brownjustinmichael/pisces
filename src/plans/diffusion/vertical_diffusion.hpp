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
#include "io/parameters.hpp"

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
			TRACE ("Operating..." << element_flags << " " << component_flags<< " " << data_in << " " << data_out);
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
			
			factory (YAML::Node i_coeff, datatype i_alpha) : alpha (i_alpha) {
				if (i_coeff.IsDefined ()) {
					coeff = i_coeff.as <datatype> ();
				} else {
					coeff = 0.0;
				}
			}
		
			virtual ~factory () {}
		
			virtual std::shared_ptr <plans::plan <datatype> > instance (plans::grid <datatype> **grids, datatype **matrices, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
				if (coeff) {
					return std::shared_ptr <plans::plan <datatype> > (new vertical_diffusion <datatype> (*grids [0], *grids [1], coeff, alpha, matrices [0], matrices [1], i_data_in, i_data_out, i_element_flags, i_component_flags));
				}
				return NULL;
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
	class background_vertical_diffusion : public implicit_plan <datatype>
	{
	public:
		background_vertical_diffusion (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, datatype i_alpha, datatype *i_diffusion, datatype *i_matrix_n, datatype *i_matrix_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) : 
		implicit_plan <datatype> (i_grid_n, i_grid_m, i_matrix_n, i_matrix_m, i_data_in, i_data_out, i_element_flags, i_component_flags),
		alpha (i_alpha),
		diffusion (i_diffusion) {
			std::stringstream debug;
			if (matrix_m) {
				for (int j = 0; j < m; ++j) {
					linalg::add_scaled (m, -diffusion [j] * alpha, grid_m.get_data (2) + j, matrix_m + j, m, m);
					// for (int i = 0; i < m; ++i) {
					// 	debug << matrix_m [i * m + j] << " ";
					// }
					// DEBUG ("MATRIX " << debug.str ());
					// debug.str ("");
				}
			} else {
				WARN ("No matrix");
			}
		}
	
		virtual ~background_vertical_diffusion () {}

		void execute () {	
			TRACE ("Operating..." << element_flags << " " << component_flags<< " " << data_in << " " << data_out);
			std::stringstream debug;
			if (*component_flags & z_solve) {
				if (1.0 - alpha != 0.0) {
					for (int j = 0; j < m; ++j) {
						linalg::matrix_matrix_multiply (1, ldn, m, diffusion [j] * (1.0 - alpha), grid_m.get_data (2) + j, data_in, 1.0, data_out + j, m, m, m);
					}
				}
			} else {
				for (int j = 0; j < m; ++j) {
					// for (int i = 0; i < ldn; ++i) {
					// 	debug << data_out [i * m + j] << " ";
					// }
					// DEBUG ("BEFORE " << debug.str ());
					// debug.str ("");
					linalg::matrix_matrix_multiply (1, ldn, m, diffusion [j], grid_m.get_data (2) + j, data_in, 1.0, data_out + j, m, m, m);
					// for (int i = 0; i < ldn; ++i) {
					// 	debug << data_out [i * m + j] << " ";
					// }
					// DEBUG ("AFTER " << debug.str ());
					// debug.str ("");
				}
			}
			TRACE ("Operation complete.");
		}

		class factory : public implicit_plan <datatype>::factory
		{
		private:
			datatype alpha;
			datatype *diffusion;
		public:
			factory (datatype i_alpha, datatype *i_diffusion) : alpha (i_alpha), diffusion (i_diffusion) {}
		
			virtual ~factory () {}
		
			virtual std::shared_ptr <plans::plan <datatype> > instance (plans::grid <datatype> **grids, datatype **matrices, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
				return std::shared_ptr <plans::plan <datatype> > (new background_vertical_diffusion <datatype> (*grids [0], *grids [1], alpha, diffusion, matrices [0], matrices [1], i_data_in, i_data_out, i_element_flags, i_component_flags));
			}
		};

		using implicit_plan <datatype>::element_flags;
		using implicit_plan <datatype>::component_flags;

	private:
		datatype alpha;
		datatype *diffusion;
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
