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

#include "../explicit_plan.hpp"

namespace plans
{
	template <class datatype>
	class horizontal_diffusion : public real_plan <datatype>
	{
	public:
		linear_diffusion (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, datatype i_coeff, datatype *i_data_source, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) : 
		real_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out, i_element_flags, i_component_flags),
		coeff (i_coeff),
		data_source (i_data_source) {
			TRACE ("Instantiating...");
		}
	
		virtual ~linear_diffusion () {}

		void execute () {	
			TRACE ("Operating...");
			for (int j = 1; j < m - 1; ++j) {
				data_out [j] += coeff * data_source [j] * ((data_in [j + 1] - data_in [j]) / (grid_m [j + 1] - grid_m [j]) - (data_in [j] - data_in [j - 1]) / (grid_m [j] - grid_m [j - 1])) / (grid_m [j + 1] - grid_m [j - 1]) * 2.0;
				data_out [j] += coeff * data_source [j] * ((data_in [m + j] - data_in [j]) / (grid_n [1] - grid_n [0]) - (data_in [j] - data_in [(n - 1) * m + j]) / (grid_n [1] - grid_n [0])) / (grid_n [1] - grid_n [0]);
				data_out [j] += coeff * (data_source [j + 1] - data_source [j - 1]) / (grid_m [j + 1] - grid_m [j - 1]) * (data_in [j + 1] - data_in [j - 1]) / (grid_m [j + 1] - grid_m [j - 1]);
				data_out [j] += coeff * (data_source [m + j] - data_source [(n - 1) * m + j]) / (grid_n [1] - grid_n [0]) * (data_in [m + j] - data_in [(n - 1) * m + j]) / (grid_n [1] - grid_n [0]);
			}
			#pragma omp parallel for
			for (int i = 1; i < n - 1; ++i) {
				for (int j = 1; j < m - 1; ++j) {
					data_out [i * m + j] += coeff * data_source [i * m + j] * ((data_in [i * m + j + 1] - data_in [i * m + j]) / (grid_m [j + 1] - grid_m [j]) - (data_in [i * m + j] - data_in [i * m + j - 1]) / (grid_m [j] - grid_m [j - 1])) / (grid_m [j + 1] - grid_m [j - 1]) * 2.0;
					data_out [i * m + j] += coeff * data_source [i * m + j] * ((data_in [(i + 1) * m + j] - data_in [i * m + j]) / (grid_n [i + 1] - grid_n [i]) - (data_in [i * m + j] - data_in [(i - 1) * m + j]) / (grid_n [i] - grid_n [i - 1])) / (grid_n [i + 1] - grid_n [i - 1]) * 2.0;
					data_out [i * m + j] += coeff * (data_source [i * m + j + 1] - data_source [i * m + j - 1]) / (grid_m [j + 1] - grid_m [j - 1]) * (data_in [i * m + j + 1] - data_in [i * m + j - 1]) / (grid_m [j + 1] - grid_m [j - 1]);
					data_out [i * m + j] += coeff * (data_source [(i + 1) * m + j] - data_source [(i - 1) * m + j]) / (grid_n [i + 1] - grid_n [i - 1]) * (data_in [(i + 1) * m + j] - data_in [(i - 1) * m + j]) / (grid_n [i + 1] - grid_n [i - 1]);
				}
			}
			for (int j = 1; j < m - 1; ++j) {
				data_out [(n - 1) * m + j] += coeff * data_source [(n - 1) * m + j] * ((data_in [(n - 1) * m + j + 1] - data_in [(n - 1) * m + j]) / (grid_m [j + 1] - grid_m [j]) - (data_in [(n - 1) * m + j] - data_in [(n - 1) * m + j - 1]) / (grid_m [j] - grid_m [j - 1])) / (grid_m [j + 1] - grid_m [j - 1]) * 2.0;
				data_out [(n - 1) * m + j] += coeff * data_source [(n - 1) * m + j] * ((data_in [j] - data_in [(n - 1) * m + j]) / (grid_n [n - 1] - grid_n [n - 2]) - (data_in [(n - 1) * m + j] - data_in [(n - 2) * m + j]) / (grid_n [n - 2] - grid_n [n - 2])) / (grid_n [n - 1] - grid_n [n - 2]);
				data_out [(n - 1) * m + j] += coeff * (data_source [(n - 1) * m + j + 1] - data_source [(n - 1) * m + j - 1]) / (grid_m [j + 1] - grid_m [j - 1]) * (data_in [(n - 1) * m + j + 1] - data_in [(n - 1) * m + j - 1]) / (grid_m [j + 1] - grid_m [j - 1]);
				data_out [(n - 1) * m + j] += coeff * (data_source [j] - data_source [(n - 2) * m + j]) / (grid_n [n - 1] - grid_n [n - 2]) * (data_in [j] - data_in [(n - 2) * m + j]) / (grid_n [n - 1] - grid_n [n - 2]);
			}
			TRACE ("Operation complete.");
		}
	
		using explicit_plan <datatype>::element_flags;
		using explicit_plan <datatype>::component_flags;

		class factory : public real_plan <datatype>::factory
		{
		private:
			datatype coeff;
		public:
			factory (datatype i_coeff) : coeff (i_coeff) {}
		
			virtual ~factory () {}
		
			virtual std::shared_ptr <plans::plan <datatype> > instance (plans::grid <datatype> **grids, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
				return std::shared_ptr <plans::plan <datatype> > (new horizontal_diffusion <datatype> (*grids [0], *grids [1], coeff, i_data_in, i_data_out, i_element_flags, i_component_flags));
			}
		};

	private:
		datatype coeff;
		datatype *data_source;
		datatype pioL2;
		using real_plan <datatype>::n;
		using real_plan <datatype>::ldn;
		using real_plan <datatype>::m;
		using real_plan <datatype>::grid_n;
		using real_plan <datatype>::data_in;
		using real_plan <datatype>::data_out;
	};
} /* plans */

#endif /* end of include guard: HORIZONTAL_DIFFUSION_HPP_E1A9847D */
