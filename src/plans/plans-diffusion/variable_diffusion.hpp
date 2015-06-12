/*!**********************************************************************
 * \file variable_diffusion.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef VARIABLE_DIFFUSION_HPP_E1A9847D
#define VARIABLE_DIFFUSION_HPP_E1A9847D

#include "linalg/utils.hpp"

#include "../real_plan.hpp"

namespace plans
{
	namespace diffusion
	{
		/*!**********************************************************************
		 * \brief A plan to enact a diffusion coefficient linearly dependent on a variable
		 ************************************************************************/
		template <class datatype>
		class linear : public real_plan <datatype>
		{
		private:
			using real_plan <datatype>::n;
			using real_plan <datatype>::ldn;
			using real_plan <datatype>::m;
			using real_plan <datatype>::grid_n;
			using real_plan <datatype>::grid_m;
			using real_plan <datatype>::data_in;
			using real_plan <datatype>::data_out;
			
			using real_plan <datatype>::element_flags;
			using real_plan <datatype>::component_flags;
			
			int bg_every, count;
			
			datatype coeff; //!< The base coefficient to multiply the source
			datatype min; //!< The minimum value of the coefficient
			datatype *data_source; //!< A pointer to the source data
			datatype pioL2; //!< The value pi / L^2 for speed
			
			std::vector <datatype> x1_vec, x2_vec, z1_vec, z2_vec, oodx_vec, oodz_vec, oodx2_vec, oodz2_vec, bg_state_vec, bg_deriv_vec;
			datatype *x1_ptr, *x2_ptr, *z1_ptr, *z2_ptr, *oodx_ptr, *oodz_ptr, *oodx2_ptr, *oodz2_ptr, *bg_state, *bg_diff, *bg_deriv;
			const datatype *pos_n, *pos_m;
			
		public:
			/*!**********************************************************************
			 * \copydoc real_plan::real_plan
			 * 
			 * \param i_coeff The base coefficient to multiply the source
			 * \param i_data_source A pointer to the source data
			 ************************************************************************/
			linear (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype i_coeff, datatype i_min, datatype *i_data_source, datatype *i_bg_diff, int i_bg_every, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) : 
			real_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out, i_element_flags, i_component_flags),
			bg_every (i_bg_every),
			coeff (i_coeff),
			min (i_min),
			data_source (i_data_source),
			bg_diff (i_bg_diff) {
				TRACE ("Instantiating...");
				pos_n = &(i_grid_n [0]);
				pos_m = &(i_grid_m [0]);
				
				x1_vec.resize (n * m);
				x1_ptr = &x1_vec [0];
				x2_vec.resize (n * m);
				x2_ptr = &x2_vec [0];
				z1_vec.resize (n * m);
				z1_ptr = &z1_vec [0];
				z2_vec.resize (n * m);
				z2_ptr = &z2_vec [0];
			
				// Calculate some derivative contributions that will be used frequently
				oodx_vec.resize (n);
				oodx_ptr = &oodx_vec [0];
				oodx2_vec.resize (n);
				oodx2_ptr = &oodx2_vec [0];
				oodz_vec.resize (m);
				oodz_ptr = &oodz_vec [0];
				oodz2_vec.resize (m);
				oodz2_ptr = &oodz2_vec [0];
				
				bg_state_vec.resize (m, 0.0);
				bg_state = &bg_state_vec [0];
				bg_deriv_vec.resize (m, 0.0);
				bg_deriv = &bg_deriv_vec [0];
				
				// Calculate 1/dx
				oodx_ptr [0] = 0.5 / (pos_n [1] - pos_n [0]);
				for (int i = 1; i < n - 1; ++i) {
					oodx_ptr [i] = 1.0 / (pos_n [i + 1] - pos_n [i - 1]);
				}
				oodx_ptr [n - 1] = 0.5 / (pos_n [n - 1] - pos_n [n - 2]);
			
				// Calculate 1/dz
				oodz_ptr [0] = 1.0 / (pos_m [1] - pos_m [0]);
				for (int i = 1; i < m - 1; ++i) {
					oodz_ptr [i] = 1.0 / (pos_m [i + 1] - pos_m [i - 1]);
				}
				oodz_ptr [m - 1] = 1.0 / (pos_m [m - 1] - pos_m [m - 2]);
				
				// Calculate 1/dx^2
				for (int i = 0; i < n - 1; ++i) {
					oodx2_ptr [i] = 1.0 / (pos_n [i + 1] - pos_n [i]);
				}
				oodx2_ptr [n - 1] = 1.0 / (pos_n [n - 1] - pos_n [n - 2]);
			
				// Calculate 1/dz&2
				for (int i = 0; i < m - 1; ++i) {
					oodz2_ptr [i] = 1.0 / (pos_m [i + 1] - pos_m [i]);
				}
				oodz2_ptr [m - 1] = 1.0 / (pos_m [m - 1] - pos_m [m - 2]);
				
				count = 0;
			}
			
			virtual ~linear () {}

			virtual int type () {
				return plan <datatype>::post;
			}
			
			/*!**********************************************************************
			 * \copydoc real_plan::execute
			 ************************************************************************/
			void execute () {
				// The full term is kappa * T * grad^2 f + kappa grad T grad f
				// div T grad f = ddx (T) ddx f + ddz (T) ddz f
				TRACE ("Operating...");
				
				// Take the second derivative of data_in in both directions
				linalg::matrix_copy (m, n - 1, data_in + m, x1_ptr);
				linalg::copy (m, data_in, x1_ptr + m * (n - 1));
				linalg::matrix_add_scaled (m, n, -1.0, data_in, x1_ptr);

				linalg::copy (m, data_in + m * (n - 1), x2_ptr);
				linalg::matrix_copy (m, n - 1, data_in, x2_ptr + m);
				linalg::matrix_add_scaled (m, n, -1.0, data_in, x2_ptr);

				linalg::scale (n, 0.0, z1_ptr + m - 1, m);
				linalg::matrix_copy (m - 1, n, data_in + 1, z1_ptr, m, m);
				linalg::matrix_add_scaled (m, n, -1.0, data_in, z1_ptr, m, m);

				linalg::matrix_copy (m, n, data_in, z2_ptr, m, m);
				linalg::matrix_add_scaled (m - 1, n, -1.0, data_in, z2_ptr + 1, m, m);

				for (int i = 0; i < n; ++i) {
					for (int j = 1; j < m - 1; ++j) {
						data_out [i * m + j] += std::max (coeff * (data_source [i * m + j] - bg_state [j]), min + 1.0e-4) * ((x1_ptr [i * m + j] * oodx2_ptr [i] - x2_ptr [i * m + j] * oodx2_ptr [(i - 1) % n]) * oodx_ptr [i] * 2.0 + (z1_ptr [i * m + j] * oodz2_ptr [j] - z2_ptr [i * m + j] * oodz2_ptr [j - 1]) * oodz_ptr [j] * 2.0);
					}
				}

				// Take the horizontal derivative of the source and data_in
				linalg::copy (m, data_source, x1_ptr + m * (n - 1));
				linalg::matrix_copy (m, n - 1, data_source + m, x1_ptr);
				linalg::matrix_add_scaled (m, n - 1, -1.0, data_source, x1_ptr + m);
				linalg::add_scaled (m, -1.0, data_source + m * (n - 1), x1_ptr);

				linalg::copy (m, data_in, x2_ptr + m * (n - 1));
				linalg::matrix_copy (m, n - 1, data_in + m, x2_ptr);
				linalg::matrix_add_scaled (m, n - 1, -1.0, data_in, x2_ptr + m);
				linalg::add_scaled (m, -1.0, data_in + m * (n - 1), x2_ptr);

				// Take the vertical derivative of the source and data_in
				linalg::scale (n, 0.0, z1_ptr + m - 1, m);
				linalg::matrix_copy (m - 1, n, data_source + 1, z1_ptr, m, m);
				linalg::matrix_add_scaled (m - 1, n, -1.0, data_source, z1_ptr + 1, m, m);

				linalg::scale (n, 0.0, z2_ptr, m);
				linalg::matrix_copy (m - 1, n, data_in + 1, z2_ptr, m, m);
				linalg::matrix_add_scaled (m - 1, n, -1.0, data_in, z2_ptr + 1, m, m);

				for (int i = 0; i < n; ++i) {
					for (int j = 1; j < m - 1; ++j) {
						data_out [i * m + j] += coeff * (x1_ptr [i * m + j] * oodx_ptr [i] * x2_ptr [i * m + j] * oodx_ptr [i] + (z1_ptr [i * m + j] - bg_deriv [j]) * oodz_ptr [j] * z2_ptr [i * m + j] * oodz_ptr [j]);
					}
				}

				if (count % bg_every == 0) {
					datatype old_bg;
					for (int j = 0; j < m; ++j) {
						old_bg = bg_state [j];
						bg_state [j] = 0.0;
						for (int i = 0; i < n; ++i) {
							bg_state [j] += data_source [i * m + j];
						}
						bg_state [j] /= n;
						bg_diff [j] += coeff * (bg_state [j] - old_bg);
					}
					bg_deriv [0] = (bg_state [1] - bg_state [0]);
					for (int j = 0; j < m; ++j) {
						bg_deriv [j] = (bg_state [j + 1] - bg_state [j - 1]);
					}
					bg_deriv [m - 1] = (bg_state [m - 1] - bg_state [m - 2]);
					*component_flags &= ~plans_setup;
				}

				count++;
				
				TRACE ("Operation complete.");
			}
			
			/*!**********************************************************************
			 * \copydoc real_plan::factory
			 ************************************************************************/
			class factory : public real_plan <datatype>::factory
			{
			private:
				datatype coeff; //!< The base coefficient for the plan to be constructed
				datatype min; //!< The value for the coefficient
				datatype *data_source; //!< The source data pointer for the plan to be constructed
				datatype *bg_diff;
				int bg_every;
				
			public:
				/*!**********************************************************************
				 * \param i_coeff The base coefficient for the plan to be constructed
				 * \param i_data_source The source data pointer for the plan to be constructed
				 ************************************************************************/
				factory (datatype i_coeff, datatype i_min, datatype *i_data_source, datatype *i_bg_diff, int i_bg_every) : coeff (i_coeff), min (i_min), data_source (i_data_source), bg_diff (i_bg_diff), bg_every (i_bg_every) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc real_plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > instance (grids::grid <datatype> **grids, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
					return std::shared_ptr <plans::plan <datatype> > (new linear <datatype> (*grids [0], *grids [1], coeff, min, data_source, bg_diff, bg_every, i_data_in, i_data_out, i_element_flags, i_component_flags));
				}
			};
		};
	} /* diffusion */
} /* plans */

#endif /* end of include guard: VARIABLE_DIFFUSION_HPP_E1A9847D */
