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
		template <class datatype>
		class variable_diffusion : public implicit_plan <datatype>
		{
		protected:
			using implicit_plan <datatype>::n;
			using implicit_plan <datatype>::ldn;
			using implicit_plan <datatype>::m;
			using implicit_plan <datatype>::grid_n;
			using implicit_plan <datatype>::grid_m;
			using implicit_plan <datatype>::data_in;
			using implicit_plan <datatype>::data_out;
			using implicit_plan <datatype>::coeff;
			using implicit_plan <datatype>::matrix_n;
			using implicit_plan <datatype>::matrix_m;
			
			using implicit_plan <datatype>::element_flags;
			using implicit_plan <datatype>::component_flags;

			datatype *data_source; //!< A pointer to the source data
			datatype *new_matrix; //!< A pointer to the diffusion matrix
			datatype *current; //!< A pointer to the current values of the diffusion
			datatype *bg_deriv; //!< A pointer to the derivative of the background source
			datatype *bg_deriv2; //!<  pointer to the second derivative of the background source
			std::vector<datatype> bg_deriv_vec; //!< A vector containing the derivative of the background diffusion
			std::vector<datatype> bg_deriv2_vec; //!< A vector containing the derivative of the background diffusion
			std::vector<datatype> new_matrix_vec; //!< A vector containing the implicit diffusion matrix
			std::vector<datatype> bg_val; //!< A vector containing the average background source
			std::vector<datatype> current_vec; //!< A vector containing the current data
			const datatype *pos_m; //!< A pointer to the vertical position, for speed
			const datatype *pos_n; //!< A pointer to the horizontal positions, for speed
			datatype *oodx; //!< An array of one over the differential of x, centered in the cell
			datatype *oodx2; //!< An array of one over the differential of x, centered at the boundary
			datatype *oodz; //!< An array of one over the differential of z, centered in the cell
			datatype *oodz2; //!< An array of one over the differential of z, centered at the boundary

		public:
			/**
			 * @param matrix_n A pointer to the horizontal matrix
			 * @param matrix_m A pointer to the vertical matrix
			 * @param i_data_source A reference to the source value that enters inside the divergence
			 * @param i_data_in A reference to the input data
			 * @param i_data_out A reference to the output data
			 * @param i_coeff The coefficient by which to multiply the results
			 */
			variable_diffusion (datatype *matrix_n, datatype *matrix_m, grids::variable <datatype> &i_data_source, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0):
			implicit_plan <datatype> (matrix_n, matrix_m, i_data_in, i_data_out, i_coeff, real_real, real_real),
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
			class factory : public implicit_plan <datatype>::factory
			{
			private:
				grids::variable <datatype> &data_source; //!< The source data pointer for the plan to be constructed
				
			public:
				/*!**********************************************************************
				 * \param i_coeff The base coefficient for the plan to be constructed
				 * \param i_data_source The source data pointer for the plan to be constructed
				 ************************************************************************/
				factory (grids::variable <datatype> &i_data_source, datatype i_coeff = 1.0) : implicit_plan <datatype>::factory (i_coeff), data_source (i_data_source) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					return std::shared_ptr <plans::plan <datatype> > (new variable_diffusion <datatype> (matrices [0], matrices [1], data_source, i_data_in, i_data_out, coeff));
				}
			};
		};

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
			
			int bg_every; //!< The number of evaluations between which to recalculate the background diffusion
			int count; //!< The current count of evaluations
			
			datatype coeff; //!< The base coefficient to multiply the source
			datatype min; //!< The minimum value of the coefficient
			datatype *data_source; //!< A pointer to the source data
			datatype pioL2; //!< The value pi / L^2 for speed
			
			std::vector <datatype> x1_vec; //!< A vector for an intermediary step
			std::vector <datatype> x2_vec; //!< A vector for an intermediary step
			std::vector <datatype> z1_vec; //!< A vector for an intermediary step
			std::vector <datatype> z2_vec; //!< A vector for an intermediary step
			std::vector <datatype> bg_state_vec; //!< A vector to contain the background state for the diffusion coefficient
			std::vector <datatype> bg_deriv_vec; //!< A vector ot contain the derivative of the background state for the diffusion coefficient

			datatype *x1_ptr; //!< A pointer to the x1 vector
			datatype *x2_ptr; //!< A pointer to the x2 vector
			datatype *z1_ptr; //!< A pointer to the z1 vector
			datatype *z2_ptr; //!< A pointer to the x2 vector
			datatype *bg_state; //!< A pointer to the bg_state vector
			datatype *bg_diff; //!< A pointer to the background diffusion coefficient
			datatype *bg_deriv; //!< A pointer to the bg_deriv vector
			const datatype *pos_m; //!< A pointer to the vertical position, for speed
			const datatype *pos_n; //!< A pointer to the horizontal positions, for speed
			datatype *oodx; //!< An array of one over the differential of x, centered in the cell
			datatype *oodx2; //!< An array of one over the differential of x, centered at the boundary
			datatype *oodz; //!< An array of one over the differential of z, centered in the cell
			datatype *oodz2; //!< An array of one over the differential of z, centered at the boundary
			
		public:
			/*!**********************************************************************
			 * \copydoc real_plan::real_plan
			 * 
			 * @param i_min The minimum value of the coefficient
			 * \param i_data_source A pointer to the source data
			 * @param i_bg_diff The background diffusion coefficient
			 * @param i_bg_every The number of evaluations between which to recalculate the background diffusion
			 ************************************************************************/
			linear (datatype i_min, grids::variable <datatype> &i_data_source, datatype *i_bg_diff, int i_bg_every, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0) : 
			real_plan <datatype> (i_data_in, i_data_out, i_coeff),
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
					component_flags &= ~plans_setup;
				}

				count++;
				
				TRACE ("Operation complete.");
			}
			
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public real_plan <datatype>::factory
			{
			protected:
				using real_plan <datatype>::factory::coeff;
				datatype min; //!< The minimum value for the coefficient
				grids::variable <datatype> &data_source; //!< The source data pointer for the plan to be constructed
				datatype *bg_diff; //!< The background diffusion coefficient
				int bg_every; //!< The number of evaluations between which to recalculate the background diffusion
				
			public:
				/*!**********************************************************************
				 * \param i_coeff The base coefficient for the plan to be constructed
				 * \param i_data_source The source data pointer for the plan to be constructed
				 * @param i_min The minimum value of the data for the plan to be constructed
				 * @param i_bg_diff The background diffusion coefficient
				 * @param i_bg_every The number of evaluations between which to recalculate the background diffusion
				 ************************************************************************/
				factory (datatype i_coeff, datatype i_min, grids::variable <datatype> &i_data_source, datatype *i_bg_diff, int i_bg_every) : real_plan <datatype>::factory (i_coeff), min (i_min), data_source (i_data_source), bg_diff (i_bg_diff), bg_every (i_bg_every) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					return std::shared_ptr <plans::plan <datatype> > (new linear <datatype> (min, data_source, bg_diff, bg_every, i_data_in, i_data_out, coeff));
				}
			};
		};
	} /* diffusion */
} /* plans */

#endif /* end of include guard: VARIABLE_DIFFUSION_HPP_E1A9847D */
