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
		/**
		 * @brief An implicit plan designed to include diffusion with an additional multiplier
		 * @details The precise form of the term is grad dot ($source) grad ($data). This usually occurs in the inclusion of density to the equations, though that is not the only case where terms like this appear. The execution of this plan does use a bit of a hack for stability: the right hand side generated is that from a full timestep ago.
		 */
		
		class variable_diffusion : public implicit_plan
		{
		protected:
			using implicit_plan::n;
			using implicit_plan::ldn;
			using implicit_plan::m;
			using implicit_plan::grid_n;
			using implicit_plan::grid_m;
			using implicit_plan::data_in;
			using implicit_plan::data_out;
			using implicit_plan::coeff;
			using implicit_plan::matrix_n;
			using implicit_plan::matrix_m;
			
			using implicit_plan::element_flags;
			using implicit_plan::component_flags;

			double *data_source; //!< A pointer to the source data
			double *new_matrix; //!< A pointer to the diffusion matrix
			double *current; //!< A pointer to the current values of the diffusion
			double *bg_deriv; //!< A pointer to the derivative of the background source
			double *bg_deriv2; //!<  pointer to the second derivative of the background source
			std::vector<double> bg_deriv_vec; //!< A vector containing the derivative of the background diffusion
			std::vector<double> bg_deriv2_vec; //!< A vector containing the derivative of the background diffusion
			std::vector<double> new_matrix_vec; //!< A vector containing the implicit diffusion matrix
			std::vector<double> bg_val; //!< A vector containing the average background source
			std::vector<double> current_vec; //!< A vector containing the current data
			const double *pos_m; //!< A pointer to the vertical position, for speed
			const double *pos_n; //!< A pointer to the horizontal positions, for speed
			double *oodx; //!< An array of one over the differential of x, centered in the cell
			double *oodx2; //!< An array of one over the differential of x, centered at the boundary
			double *oodz; //!< An array of one over the differential of z, centered in the cell
			double *oodz2; //!< An array of one over the differential of z, centered at the boundary

		public:
			/**
			 * @param matrix_n A pointer to the horizontal matrix
			 * @param matrix_m A pointer to the vertical matrix
			 * @param i_data_source A reference to the source value that enters inside the divergence
			 * @param i_data_in A reference to the input data
			 * @param i_data_out A reference to the output data
			 * @param i_coeff The coefficient by which to multiply the results
			 */
			variable_diffusion (double *matrix_n, double *matrix_m, grids::variable &i_data_source, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0):
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
			
			virtual ~variable_diffusion () {}

			/**
			 * @copydoc plan::setup
			 */
			void setup () {
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

			/**
			 * @copydoc plan::execute
			 */
			virtual void execute () {
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

			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public implicit_plan::factory
			{
			private:
				grids::variable &data_source; //!< The source data pointer for the plan to be constructed
				
			public:
				/*!**********************************************************************
				 * \param i_coeff The base coefficient for the plan to be constructed
				 * \param i_data_source The source data pointer for the plan to be constructed
				 ************************************************************************/
				factory (grids::variable &i_data_source, double i_coeff = 1.0) : implicit_plan::factory (i_coeff), data_source (i_data_source) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					return std::shared_ptr <plans::plan > (new variable_diffusion (matrices [0], matrices [1], data_source, i_data_in, i_data_out, coeff));
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to enact a diffusion coefficient linearly dependent on a variable
		 ************************************************************************/
		
		class linear : public real_plan
		{
		protected:
			using real_plan::n;
			using real_plan::ldn;
			using real_plan::m;
			using real_plan::grid_n;
			using real_plan::grid_m;
			using real_plan::data_in;
			using real_plan::data_out;
			
			using real_plan::element_flags;
			using real_plan::component_flags;
			
			int bg_every; //!< The number of evaluations between which to recalculate the background diffusion
			int count; //!< The current count of evaluations
			
			double coeff; //!< The base coefficient to multiply the source
			double min; //!< The minimum value of the coefficient
			double *data_source; //!< A pointer to the source data
			double pioL2; //!< The value pi / L^2 for speed
			
			std::vector <double> x1_vec; //!< A vector for an intermediary step
			std::vector <double> x2_vec; //!< A vector for an intermediary step
			std::vector <double> z1_vec; //!< A vector for an intermediary step
			std::vector <double> z2_vec; //!< A vector for an intermediary step
			std::vector <double> bg_state_vec; //!< A vector to contain the background state for the diffusion coefficient
			std::vector <double> bg_deriv_vec; //!< A vector ot contain the derivative of the background state for the diffusion coefficient

			double *x1_ptr; //!< A pointer to the x1 vector
			double *x2_ptr; //!< A pointer to the x2 vector
			double *z1_ptr; //!< A pointer to the z1 vector
			double *z2_ptr; //!< A pointer to the x2 vector
			double *bg_state; //!< A pointer to the bg_state vector
			double *bg_diff; //!< A pointer to the background diffusion coefficient
			double *bg_deriv; //!< A pointer to the bg_deriv vector
			const double *pos_m; //!< A pointer to the vertical position, for speed
			const double *pos_n; //!< A pointer to the horizontal positions, for speed
			double *oodx; //!< An array of one over the differential of x, centered in the cell
			double *oodx2; //!< An array of one over the differential of x, centered at the boundary
			double *oodz; //!< An array of one over the differential of z, centered in the cell
			double *oodz2; //!< An array of one over the differential of z, centered at the boundary
			
		public:
			/*!**********************************************************************
			 * \copydoc real_plan::real_plan
			 * 
			 * @param i_min The minimum value of the coefficient
			 * \param i_data_source A pointer to the source data
			 * @param i_bg_diff The background diffusion coefficient
			 * @param i_bg_every The number of evaluations between which to recalculate the background diffusion
			 ************************************************************************/
			linear (double i_min, grids::variable &i_data_source, double *i_bg_diff, int i_bg_every, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0) : 
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
			
			virtual ~linear () {}

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

			void reset_source (double *new_source) {
				data_source = new_source;
			}
			
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public real_plan::factory
			{
			protected:
				using real_plan::factory::coeff;
				double min; //!< The minimum value for the coefficient
				grids::variable &data_source; //!< The source data pointer for the plan to be constructed
				double *bg_diff; //!< The background diffusion coefficient
				int bg_every; //!< The number of evaluations between which to recalculate the background diffusion
				
			public:
				/*!**********************************************************************
				 * \param i_coeff The base coefficient for the plan to be constructed
				 * \param i_data_source The source data pointer for the plan to be constructed
				 * @param i_min The minimum value of the data for the plan to be constructed
				 * @param i_bg_diff The background diffusion coefficient
				 * @param i_bg_every The number of evaluations between which to recalculate the background diffusion
				 ************************************************************************/
				factory (double i_coeff, double i_min, grids::variable &i_data_source, double *i_bg_diff, int i_bg_every) : real_plan::factory (i_coeff), min (i_min), data_source (i_data_source), bg_diff (i_bg_diff), bg_every (i_bg_every) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					return std::shared_ptr <plans::plan > (new linear (min, data_source, bg_diff, bg_every, i_data_in, i_data_out, coeff));
				}
			};
		};

		class tanh : public linear
		{
		protected:
			using linear::n;
			using linear::ldn;
			using linear::m;
			using linear::grid_n;
			using linear::grid_m;
			using linear::data_in;
			using linear::data_out;
			
			using linear::element_flags;
			using linear::component_flags;
			
			double source_length;
			double source_zero;

			std::vector <double> src_vec;
			double *src_ptr;
			double *data_orig; //!< A pointer to the source data

		public:
			/*!**********************************************************************
			 * \copydoc real_plan::real_plan
			 * 
			 * @param i_min The minimum value of the coefficient
			 * \param i_data_source A pointer to the source data
			 * @param i_bg_diff The background diffusion coefficient
			 * @param i_bg_every The number of evaluations between which to recalculate the background diffusion
			 ************************************************************************/
			tanh (double i_source_length, double i_source_zero, double i_min, grids::variable &i_data_source, double *i_bg_diff, int i_bg_every, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0) : 
			linear (i_min, i_data_source, i_bg_diff, i_bg_every, i_data_in, i_data_out, i_coeff),
			source_length (i_source_length),
			source_zero (i_source_zero),
			data_orig (i_data_source.ptr ()) {
				TRACE ("Instantiating...");
				src_vec.resize (n * m);
				src_ptr = &src_vec [0];
				linear::reset_source (src_ptr);
			}
			
			virtual ~tanh () {}

			/*!**********************************************************************
			 * \copydoc real_plan::execute
			 ************************************************************************/
			void execute () {
				// The full term is kappa * T * grad^2 f + kappa grad T grad f
				// div T grad f = ddx (T) ddx f + ddz (T) ddz f
				TRACE ("Operating...");
				
				// Calculate the source value
				int g;
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j < m; ++j) {
						g = i * m + j;
						if (data_source [g] > source_zero)
						{
							src_ptr [g] = std::tanh ((data_orig [g] - source_zero) / source_length);
						} else {
							src_ptr [g] = std::tanh ((data_orig [g] - source_zero) / source_length);
						}
					}
				}

				linear::execute ();
				TRACE ("Operation complete.");
			}
			
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public real_plan::factory
			{
			protected:
				using real_plan::factory::coeff;
				double min; //!< The minimum value for the coefficient
				double source_length;
				double source_zero;
				grids::variable &data_source; //!< The source data pointer for the plan to be constructed
				double *bg_diff; //!< The background diffusion coefficient
				int bg_every; //!< The number of evaluations between which to recalculate the background diffusion
				
			public:
				/*!**********************************************************************
				 * \param i_coeff The base coefficient for the plan to be constructed
				 * \param i_data_source The source data pointer for the plan to be constructed
				 * @param i_min The minimum value of the data for the plan to be constructed
				 * @param i_bg_diff The background diffusion coefficient
				 * @param i_bg_every The number of evaluations between which to recalculate the background diffusion
				 ************************************************************************/
				factory (double i_coeff, double i_source_length, double i_source_zero, double i_min, grids::variable &i_data_source, double *i_bg_diff, int i_bg_every) : real_plan::factory (i_coeff), min (i_min), source_length (i_source_length), source_zero (i_source_zero), data_source (i_data_source), bg_diff (i_bg_diff), bg_every (i_bg_every) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					return std::shared_ptr <plans::plan > (new tanh (source_length, source_zero, min, data_source, bg_diff, bg_every, i_data_in, i_data_out, coeff));
				}
			};
		};
	} /* diffusion */
} /* plans */

#endif /* end of include guard: VARIABLE_DIFFUSION_HPP_E1A9847D */
