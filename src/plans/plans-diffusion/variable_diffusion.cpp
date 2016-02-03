/*!**********************************************************************
 * \file variable_diffusion.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "variable_diffusion.hpp"

namespace plans
{
	namespace diffusion
	{
		/**
		 * @brief An implicit plan designed to include diffusion with an additional multiplier
		 * @details The precise form of the term is grad dot ($source) grad ($data). This usually occurs in the inclusion of density to the equations, though that is not the only case where terms like this appear. The execution of this plan does use a bit of a hack for stability: the right hand side generated is that from a full timestep ago.
		 */
		
		variable_diffusion::variable_diffusion (double *matrix_n, double *matrix_m, grids::variable &i_data_source, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff):
		implicit_plan (matrix_n, matrix_m, i_data_in, i_data_out, i_coeff, real_real, real_real),
		data_source (i_data_source.ptr (real_real)) {
			pos_n = &grid_n [0];
			pos_m = &grid_m [0];
			bg_deriv_vec.resize (m, 0.0);
			bg_deriv = &bg_deriv_vec [0];
			bg_deriv2_vec.resize (m, 0.0);
			bg_deriv2 = &bg_deriv2_vec [0];
			bg_val.resize (m, 0.0);
			new_matrix_vec.resize (m * m, 0.0);
			new_matrix = &new_matrix_vec [0];
			current_vec.resize (m * ldn, 0.0);
			current = &current_vec [0];

			oodx = grid_n.get_ood ();
			oodx2 = grid_n.get_ood2 ();
			oodz = grid_m.get_ood ();
			oodz2 = grid_m.get_ood2 ();
		}
		
		void variable_diffusion::setup () {
			TRACE ("Setting up with coefficient " << coeff);
			if (matrix_m) {
				linalg::add_scaled (m, coeff * log (data_source [1] / data_source [0]) / (pos_m [1] - pos_m [0]), grid_m.get_data (1), new_matrix, m, m);
				for (int j = 1; j < m - 1; ++j) {
					bg_deriv [j] = 0.0;
					bg_deriv2 [j] = 0.0;
					for (int i = 0; i < n; ++i)
					{
						bg_deriv [j] += log (data_source [i * m + j + 1] / data_source [i * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]);
						bg_deriv2 [j] += log (data_source [i * m + j + 1] / data_source [i * m + j]) / (pos_m [j + 1] - pos_m [j]) - log (data_source [i * m + j] / data_source [i * m + j - 1]) / (pos_m [j] - pos_m [j - 1]);
					}
					bg_deriv [j] /= n;
					bg_deriv2 [j] /= n;
				}

				for (int j = 1; j < m - 1; ++j) {
					linalg::add_scaled (m, coeff * bg_deriv [j], grid_m.get_data (1) + j, new_matrix + j, m, m);

					linalg::add_scaled (m, coeff * bg_deriv2 [j] / (pos_m [j + 1] - pos_m [j - 1]) * 2., grid_m.get_data (0) + j, new_matrix + j, m, m);
				}
				linalg::add_scaled (m, coeff * log (data_source [m - 1] / data_source [m - 2]), grid_m.get_data (1) + m - 1, new_matrix + m - 1, m, m);
				linalg::add_scaled (m * m, -1.0, new_matrix, matrix_m);
			} else {
				WARN ("No matrix");
			}
		}

		void variable_diffusion::execute () {
			linalg::add_scaled (m * ldn, 1.0, current, data_out);
			linalg::scale (m * ldn, 0.0, current);
			#pragma omp parallel for
			for (int i = 0; i < n; ++i)
			{
				int m1 = 0, p1 = 0, g = 0;
				p1 = (i + 1) % n;
				m1 = (i - 1 + n) % n;
				for (int j = 1; j < m - 1; ++j)
				{
					g = i * m + j;
					current [g] += coeff * (log (data_source [g + 1] / data_source [g - 1]) * oodz2 [j] - bg_deriv [j]) * (data_in [g + 1] - data_in [g - 1]) * oodz2 [j];

					current [g] += coeff * (log (data_source [p1 * m + j] / data_source [m1 * m + j]) * oodx2 [i]) * (data_in [p1 * m + j] - data_in [m1 * m + j]) * oodx2 [i];
				}
			}
		}

		std::shared_ptr <plans::plan > variable_diffusion::factory::_instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
			return std::shared_ptr <plans::plan > (new variable_diffusion (matrices [0], matrices [1], data_source, i_data_in, i_data_out, coeff));
		}

		linear::linear (double i_min, grids::variable &i_data_source, double *i_bg_diff, int i_bg_every, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff) : 
		real_plan (i_data_in, i_data_out, i_coeff),
		bg_every (i_bg_every),
		coeff (i_coeff),
		min (i_min),
		data_source (i_data_source.ptr ()),
		bg_diff (i_bg_diff) {
			TRACE ("Instantiating...");
			pos_n = &(grid_n [0]);
			pos_m = &(grid_m [0]);
			
			x1_vec.resize (n * m);
			x1_ptr = &x1_vec [0];
			x2_vec.resize (n * m);
			x2_ptr = &x2_vec [0];
			z1_vec.resize (n * m);
			z1_ptr = &z1_vec [0];
			z2_vec.resize (n * m);
			z2_ptr = &z2_vec [0];
		
			// Calculate some derivative contributions that will be used frequently
			oodx = grid_n.get_ood ();
			oodx2 = grid_n.get_ood2 ();
			oodz = grid_m.get_ood ();
			oodz2 = grid_m.get_ood2 ();
			
			bg_state_vec.resize (m, 0.0);
			bg_state = &bg_state_vec [0];
			bg_deriv_vec.resize (m, 0.0);
			bg_deriv = &bg_deriv_vec [0];
			
			count = 0;
		}
			
		void linear::execute () {
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
					data_out [i * m + j] += std::max (coeff * (data_source [i * m + j] - bg_state [j]), min + 1.0e-4) * ((x1_ptr [i * m + j] * oodx2 [i] - x2_ptr [i * m + j] * oodx2 [(i - 1) % n]) * oodx [i] * 2.0 + (z1_ptr [i * m + j] * oodz2 [j] - z2_ptr [i * m + j] * oodz2 [j - 1]) * oodz [j] * 2.0);
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
					data_out [i * m + j] += coeff * (x1_ptr [i * m + j] * oodx [i] * x2_ptr [i * m + j] * oodx [i] + (z1_ptr [i * m + j] - bg_deriv [j]) * oodz [j] * z2_ptr [i * m + j] * oodz [j]);
				}
			}

			if (count % bg_every == 0) {
				double old_bg;
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
				for (int j = 1; j < m - 1; ++j) {
					bg_deriv [j] = (bg_state [j + 1] - bg_state [j - 1]);
				}
				bg_deriv [m - 1] = (bg_state [m - 1] - bg_state [m - 2]);
				component_flags &= ~plans_setup;
			}

			count++;
			
			TRACE ("Operation complete.");
		}

		void linear::reset_source (double *new_source) {
			data_source = new_source;
		}
		
		std::shared_ptr <plans::plan > linear::factory::_instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
			return std::shared_ptr <plans::plan > (new linear (min, data_source, bg_diff, bg_every, i_data_in, i_data_out, coeff));
		}

		tanh::tanh (double i_source_length, double i_source_zero, double i_min, grids::variable &i_data_source, double *i_bg_diff, int i_bg_every, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff) : 
		linear (i_min, i_data_source, i_bg_diff, i_bg_every, i_data_in, i_data_out, i_coeff),
		source_length (i_source_length),
		source_zero (i_source_zero),
		data_orig (i_data_source.ptr ()) {
			TRACE ("Instantiating...");
			src_vec.resize (n * m);
			src_ptr = &src_vec [0];
			linear::reset_source (src_ptr);
		}
			
		void tanh::execute () {
			// The full term is kappa * T * grad^2 f + kappa grad T grad f
			// div T grad f = ddx (T) ddx f + ddz (T) ddz f
			TRACE ("Operating...");
			
			// Calculate the source value
			int g;
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < m; ++j) {
					g = i * m + j;
					src_ptr [g] = std::tanh ((data_orig [g] - source_zero) / source_length);
				}
			}

			linear::execute ();
			TRACE ("Operation complete.");
		}

		std::shared_ptr <plans::plan > tanh::factory::_instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
			return std::shared_ptr <plans::plan > (new tanh (source_length, source_zero, min, data_source, bg_diff, bg_every, i_data_in, i_data_out, coeff));
		}

		ra_tanh::ra_tanh (double i_source_length, double i_source_zero, double chi, double stiffness, double i_min, grids::variable &i_data_source, double *i_bg_diff, int i_bg_every, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff) : 
		linear (i_min, i_data_source, i_bg_diff, i_bg_every, i_data_in, i_data_out, i_coeff),
		source_length (i_source_length),
		source_zero (i_source_zero),
		data_orig (i_data_source.ptr ()) {
			TRACE ("Instantiating...");
			src_vec.resize (n * m);
			src_ptr = &src_vec [0];
			linear::reset_source (src_ptr);
			a = chi - 1.;
			b = 1. + 2. * (stiffness - 2.) * chi + chi * chi;
			c = 2. * (1. + stiffness) * chi;
			d = chi * (1. - stiffness);
			e = chi * (1 + stiffness);
		}
			
		void ra_tanh::execute () {
			// The full term is kappa * T * grad^2 f + kappa grad T grad f
			// div T grad f = ddx (T) ddx f + ddz (T) ddz f
			TRACE ("Operating...");
			
			// Calculate the source value
			int g;
			double arg;
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < m; ++j) {
					g = i * m + j;
					arg = -(data_orig [g] - source_zero) / source_length;
					src_ptr [g] = (a - sqrt(b - c * std::tanh(arg))) / (d + e * std::tanh(arg));
					DEBUG (b - c * std::tanh(arg) << " " << src_ptr [g]);
				}
			}

			linear::execute ();
			TRACE ("Operation complete.");
		}
			
		std::shared_ptr <plans::plan > ra_tanh::factory::_instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
			return std::shared_ptr <plans::plan > (new ra_tanh (source_length, source_zero, phi, stiffness, min, data_source, bg_diff, bg_every, i_data_in, i_data_out, coeff));
		}
	} /* diffusion */
} /* plans */
