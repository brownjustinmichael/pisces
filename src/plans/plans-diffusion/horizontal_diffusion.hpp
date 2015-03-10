/*!**********************************************************************
 * \file horizontal_diffusion.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef HORIZONTAL_DIFFUSION_HPP_E1A9847D
#define HORIZONTAL_DIFFUSION_HPP_E1A9847D

#include "linalg/utils.hpp"

#include "../implicit_plan.hpp"
#include "io/parameters.hpp"

namespace plans
{
	template <class datatype>
	class horizontal_diffusion : public implicit_plan <datatype>
	{
	public:
		horizontal_diffusion (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype i_coeff, datatype i_alpha, datatype *i_matrix_n, datatype *i_matrix_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) : 
		implicit_plan <datatype> (i_grid_n, i_grid_m, i_matrix_n, i_matrix_m, i_data_in, i_data_out, i_element_flags, i_component_flags),
		coeff (i_coeff),
		alpha (i_alpha) {
			TRACE ("Instantiating...");
			pioL2 = 4.0 * (std::acos (-1.0) * std::acos (-1.0) / (grid_n [n - 1] - grid_n [0]) / (grid_n [n - 1] - grid_n [0]));
			if (matrix_n) {
				for (int j = 0; j < m; ++j) {
					matrix_n [j] = matrix_n [m + j] = 0.0;
					for (int i = 2; i < ldn; ++i) {
						matrix_n [i * m + j] = coeff * alpha * pioL2 * (datatype) ((i / 2) * (i / 2));
					}
				}
			} else {
				WARN ("No matrix");
			}
		}
	
		virtual ~horizontal_diffusion () {}

		void execute () {	
			TRACE ("Operating..." << element_flags);
			if (*component_flags & x_solve) {
				if (1.0 - alpha != 0.0) {
					// #pragma omp parallel for
					for (int i = 2; i < ldn; ++i) {
						linalg::add_scaled (m, -coeff * (1.0 - alpha) * pioL2 * (i / 2) * (i / 2), data_in + i * m, data_out + i * m);
					}
				}
			} else {
				// #pragma omp parallel for
				for (int i = 2; i < ldn; ++i) {
					linalg::add_scaled (m, -coeff * pioL2 * (i / 2) * (i / 2), data_in + i * m, data_out + i * m);
				}
			}
			TRACE ("Operation complete.");
		}
	
		using implicit_plan <datatype>::element_flags;
		using implicit_plan <datatype>::component_flags;

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
		
			virtual std::shared_ptr <plans::plan <datatype> > instance (grids::grid <datatype> **grids, datatype **matrices, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
				if (coeff) {
					return std::shared_ptr <plans::plan <datatype> > (new horizontal_diffusion <datatype> (*grids [0], *grids [1], coeff, alpha, matrices [0], matrices [1], i_data_in, i_data_out, i_element_flags, i_component_flags));
				}
				return NULL;
			}
		};

	private:
		datatype coeff;
		datatype alpha;
		datatype pioL2;
		using implicit_plan <datatype>::n;
		using implicit_plan <datatype>::ldn;
		using implicit_plan <datatype>::m;
		using implicit_plan <datatype>::grid_n;
		using implicit_plan <datatype>::data_in;
		using implicit_plan <datatype>::data_out;
		using implicit_plan <datatype>::matrix_n;
	};
	
	template <class datatype>
	class background_horizontal_diffusion : public implicit_plan <datatype>
	{
	public:
		background_horizontal_diffusion (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype i_alpha, datatype *i_diffusion, datatype *i_matrix_n, datatype *i_matrix_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) : 
		implicit_plan <datatype> (i_grid_n, i_grid_m, i_matrix_n, i_matrix_m, i_data_in, i_data_out, i_element_flags, i_component_flags),
		alpha (i_alpha),
		diffusion (i_diffusion) {
			TRACE ("Instantiating...");
			pioL2 = 4.0 * (std::acos (-1.0) * std::acos (-1.0) / (grid_n [n - 1] - grid_n [0]) / (grid_n [n - 1] - grid_n [0]));
			if (matrix_n) {
				for (int j = 0; j < m; ++j) {
					matrix_n [j] = matrix_n [m + j] = 0.0;
					for (int i = 2; i < ldn; ++i) {
						matrix_n [i * m + j] = diffusion [j] * alpha * pioL2 * (datatype) ((i / 2) * (i / 2));
					}
				}
			} else {
				WARN ("No matrix");
			}
		}
	
		virtual ~background_horizontal_diffusion () {}

		void execute () {	
			TRACE ("Operating..." << element_flags);
			if (*component_flags & x_solve) {
				if (1.0 - alpha != 0.0) {
					// #pragma omp parallel for
					for (int j = 0; j < m; ++j) {
						for (int i = 2; i < ldn; ++i) {
							data_out [i * m + j] -= diffusion [j] * (1.0 - alpha) * pioL2 * (i / 2) * (i / 2) * data_in [i * m + j];
						}
					}
				}
			} else {
				// #pragma omp parallel for
				for (int j = 0; j < m; ++j) {
					for (int i = 2; i < ldn; ++i) {
						data_out [i * m + j] -= diffusion [j] * pioL2 * (i / 2) * (i / 2) * data_in [i * m + j];
					}
				}
			}
			TRACE ("Operation complete.");
		}
	
		using implicit_plan <datatype>::element_flags;
		using implicit_plan <datatype>::component_flags;

		class factory : public implicit_plan <datatype>::factory
		{
		private:
			datatype alpha, *diffusion;
		public:
			factory (datatype i_alpha, datatype *i_diffusion) : alpha (i_alpha), diffusion (i_diffusion) {}
		
			virtual ~factory () {}
		
			virtual std::shared_ptr <plans::plan <datatype> > instance (grids::grid <datatype> **grids, datatype **matrices, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
				return std::shared_ptr <plans::plan <datatype> > (new background_horizontal_diffusion <datatype> (*grids [0], *grids [1], alpha, diffusion, matrices [0], matrices [1], i_data_in, i_data_out, i_element_flags, i_component_flags));
			}
		};

	private:
		datatype alpha;
		datatype *diffusion;
		datatype pioL2;
		using implicit_plan <datatype>::n;
		using implicit_plan <datatype>::ldn;
		using implicit_plan <datatype>::m;
		using implicit_plan <datatype>::grid_n;
		using implicit_plan <datatype>::data_in;
		using implicit_plan <datatype>::data_out;
		using implicit_plan <datatype>::matrix_n;
	};
} /* plans */

#endif /* end of include guard: HORIZONTAL_DIFFUSION_HPP_E1A9847D */
