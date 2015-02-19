/*!**********************************************************************
 * \file advection_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces/src/two_d
 * 
 * Created by Justin Brown on 2013-12-16.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ADVECTION_TWO_D_HPP_GGR0NN1Q
#define ADVECTION_TWO_D_HPP_GGR0NN1Q

#include "../real_plan.hpp"
#include "io/parameters.hpp"

namespace plans
{
	void advection_term (int i, int j, int n, int m, double coeff, const double *pos_n, const double *pos_m, const double *vel_n, const double *vel_m, const double *data_in, double *data_out) {
		data_out [i * m + j] += coeff * (vel_n [i * m + j] * (data_in [((i + 1) % n) * m + j] - data_in [((i - 1) % n) * m + j]) / (pos_n [(i + 1) % n] - pos_n [(i - 1) % n]) + vel_m [i * m + j] * (data_in [i * m + j + 1] - data_in [i * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]));
	}

	void advection_term_quadratic (int i, int j, int n, int m, double coeff, const double *pos_n, const double *pos_m, const double *vel_n, const double *vel_m, const double *data_in, double *data_out) {
		double denom = (pos_n [(i + 1) % n] - pos_n [i]) * (pos_n [(i + 1) % n] - pos_n [(i - 1) % n]) * (pos_n [i] - pos_n [(i - 1) % n]);
		double A = (pos_n [(i - 1) % n] * (data_in [i * m + j] - data_in [((i + 1) % n) * m + j]) + pos_n [i] * (data_in [((i + 1) % n) * m + j] - data_in [((i - 1) % n) * m + j]) + pos_n [(i + 1) % n] * (data_in [((i - 1) % n) * m + j] - data_in [i * m + j])) / denom;
		double B = (pos_n [(i - 1) % n] * pos_n [(i - 1) % n] * (data_in [((i + 1) % n) * m + j] - data_in [i * m + j]) + pos_n [i] * pos_n [i] * (data_in [((i - 1) % n) * m + j] - data_in [((i + 1) % n) * m + j]) + pos_n [(i + 1) % n] * pos_n [(i + 1) % n] * (data_in [i * m + j] - data_in [((i - 1) % n) * m + j])) / denom;
		data_out [i * m + j] += coeff * vel_n [i * m + j] * (2.0 * A * pos_n [i] + B);
		
		denom = (pos_m [j + 1] - pos_m [j]) * (pos_m [j + 1] - pos_m [j - 1]) * (pos_m [j] - pos_m [j - 1]);
		A = (pos_m [j - 1] * (data_in [i * m + j] - data_in [i * m + j + 1]) + pos_m [j] * (data_in [i * m + j + 1] - data_in [i * m + j - 1]) + pos_m [j + 1] * (data_in [i * m + j - 1] - data_in [i * m + j])) / denom;
		B = (pos_m [j - 1] * pos_m [j - 1] * (data_in [i * m + j + 1] - data_in [i * m + j]) + pos_m [j] * pos_m [j] * (data_in [i * m + j - 1] - data_in [i * m + j + 1]) + pos_m [j + 1] * pos_m [j + 1] * (data_in [i * m + j] - data_in [i * m + j - 1])) / denom;
		data_out [i * m + j] += coeff * vel_m [i * m + j] * (2.0 * A * pos_m [j] + B);
	}

	template <class datatype>
	class advection : public real_plan <datatype>
	{
	public:
		advection (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, datatype i_coeff, datatype* i_vel_n, datatype *i_vel_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) :
		real_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out, i_element_flags, i_component_flags),
		coeff (i_coeff),
		vel_n (i_vel_n),
		vel_m (i_vel_m),
		pos_n (&(grid_n [0])),
		pos_m (&(grid_m [0])) {
			TRACE ("Adding advection...");
			x_vec.resize (n * m);
			x_ptr = &x_vec [0];
			z_vec.resize (n * m);
			z_ptr = &z_vec [0];
			oodx_vec.resize (n);
			oodx_ptr = &oodx_vec [0];
			oodz_vec.resize (m);
			oodz_ptr = &oodz_vec [0];
			
			oodx_ptr [0] = 1.0 / (pos_n [1] - pos_n [0]);
			for (int i = 1; i < n - 1; ++i) {
				oodx_ptr [i] = 1.0 / (pos_n [i + 1] - pos_n [i - 1]);
			}
			oodx_ptr [n - 1] = 1.0 / (pos_n [n - 1] - pos_n [n - 2]);
			
			oodz_ptr [0] = 1.0 / (pos_m [1] - pos_m [0]);
			for (int i = 1; i < m - 1; ++i) {
				oodz_ptr [i] = 1.0 / (pos_m [i + 1] - pos_m [i - 1]);
			}
			oodz_ptr [m - 1] = 1.0 / (pos_m [m - 1] - pos_m [m - 2]);
		}
	
		virtual ~advection () {}
	
		virtual void execute () {
			linalg::matrix_copy (m, n - 1, data_in, x_ptr + m);
			linalg::matrix_add_scaled (m, n - 1, -1.0, data_in + m, x_ptr);
			linalg::add_scaled (m, -1.0, data_in, x_ptr);
			linalg::add_scaled (m, 1.0, data_in + m * (n - 1), x_ptr + m * (n - 1));

			linalg::matrix_copy (m - 1, n, data_in, z_ptr + 1);
			linalg::matrix_add_scaled (m - 1, n, -1.0, data_in + 1, z_ptr);
			linalg::add_scaled (n, -1.0, data_in + m - 1, z_ptr, m, m);
			linalg::add_scaled (n, 1.0, data_in, z_ptr + m - 1, m, m);
			
			#pragma omp parallel for
			for (int j = 0; j < m; ++j) {
				x_ptr [j] = 0.0;
				x_ptr [(n - 1) * m + j] = 0.0;
				for (int i = 1; i < n - 1; ++i) {
					x_ptr [i * m + j] = (vel_n [i * m + j] * x_ptr [i * m + j] * oodx_ptr [i] + vel_m [i * m + j] * z_ptr [i * m + j] * oodz_ptr [j]);
				}
			}

			linalg::matrix_add_scaled (m, n, coeff, x_ptr, data_out);
		}
	
		class factory : public real_plan <datatype>::factory
		{
		private:
			datatype coeff;
			datatype *vel_n, *vel_m;
		public:
			factory (datatype i_coeff, datatype* i_vel_n, datatype* i_vel_m) : coeff (i_coeff), vel_n (i_vel_n), vel_m (i_vel_m) {}

			factory (YAML::Node i_coeff, datatype* i_vel_n, datatype* i_vel_m) : vel_n (i_vel_n), vel_m (i_vel_m) {
				if (i_coeff.IsDefined ()) {
					coeff = i_coeff.as <datatype> ();
				} else {
					coeff = 0.0;
				}
			}
		
			virtual ~factory () {}
		
			virtual std::shared_ptr <plans::plan <datatype> > instance (plans::grid <datatype> **grids, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
				if (coeff) {
					return std::shared_ptr <plans::plan <datatype> > (new advection <datatype> (*grids [0], *grids [1], coeff, vel_n, vel_m, i_data_in, i_data_out, i_element_flags, i_component_flags));
				}
				return NULL;
			}
		};

	private:
		datatype coeff;
		datatype *vel_n, *vel_m;
		const datatype *pos_n, *pos_m;
		std::vector <datatype> x_vec, z_vec, oodx_vec, oodz_vec;
		datatype *x_ptr, *z_ptr, *oodx_ptr, *oodz_ptr;
		using real_plan <datatype>::n;
		using real_plan <datatype>::m;
		using real_plan <datatype>::grid_n;
		using real_plan <datatype>::grid_m;
		using real_plan <datatype>::data_in;
		using real_plan <datatype>::data_out;
	};

	template <class datatype>
	class stream_advection : public real_plan <datatype>
	{
	public:
		stream_advection (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, datatype i_coeff, datatype* i_stream, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) :
		real_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out, i_element_flags, i_component_flags),
		coeff (-i_coeff),
		stream (i_stream),
		pos_n (&(grid_n [0])),
		pos_m (&(grid_m [0])) {}
	
		virtual ~stream_advection () {}
	
		virtual void execute () {
			for (int j = 1; j < m - 1; ++j) {
				data_out [j] += coeff * ((stream [m + j] - stream [j]) * (data_in [j + 1] - data_in [j - 1]) - (stream [j + 1] - stream [j - 1]) * (data_in [m + j] - data_in [j])) / (pos_n [1] - pos_n [0]) / (pos_m [j + 1] - pos_m [j - 1]);
				for (int i = 1; i < n - 1; ++i) {
					data_out [i * m + j] += coeff * ((stream [(i + 1) * m + j] - stream [(i - 1) * m + j]) * (data_in [i * m + j + 1] - data_in [i * m + j - 1]) - (stream [i * m + j + 1] - stream [i * m + j - 1]) * (data_in [(i + 1) * m + j] - data_in [(i - 1) * m + j])) / (pos_m [j + 1] - pos_m [j - 1]) / (pos_n [i + 1] - pos_n [i - 1]);
				}
				data_out [(n - 1) * m + j] += coeff * ((stream [(n - 1) * m + j] - stream [(n - 2) * m + j]) * (data_in [(n - 1) * m + j + 1] - data_in [(n - 1) * m + j - 1]) - (stream [(n - 1) * m + j + 1] - stream [(n - 1) * m + j - 1]) * (data_in [(n - 1) * m + j] - data_in [(n - 2) * m + j])) / (pos_n [n - 1] - pos_n [n - 2]) / (pos_m [j + 1] - pos_m [j - 1]);
			}
		}

	private:
		datatype coeff;
		datatype *stream;
		const datatype *pos_n, *pos_m;
		using real_plan <datatype>::n;
		using real_plan <datatype>::m;
		using real_plan <datatype>::grid_n;
		using real_plan <datatype>::grid_m;
		using real_plan <datatype>::data_in;
		using real_plan <datatype>::data_out;
	};
} /* plans */

#endif /* end of include guard: ADVECTION_TWO_D_HPP_GGR0NN1Q */
