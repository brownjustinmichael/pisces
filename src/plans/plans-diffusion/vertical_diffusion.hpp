/*!**********************************************************************
 * \file vertical_diffusion.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DIFFUSION_TWO_D_HPP_YHECX9VS
#define DIFFUSION_TWO_D_HPP_YHECX9VS

#include <ctime>
#include <chrono>

#include "../implicit_plan.hpp"
#include "../explicit_plan.hpp"
#include "../real_plan.hpp"
#include "io/parameters.hpp"

/*!**********************************************************************
 * \namespace plans::diffusion
 * 
 * \brief A namespace containing the diffusion plan extension
 ************************************************************************/
namespace plans
{
	namespace diffusion
	{
		/*!**********************************************************************
		 * \brief An implicit plan for uniform diffusion in the z direction
		 ************************************************************************/
		class vertical : public implicit_plan
		{
		protected:
			using implicit_plan::coeff;
			using implicit_plan::n;
			using implicit_plan::ldn;
			using implicit_plan::m;
			using implicit_plan::dims;
			using implicit_plan::data_in;
			using implicit_plan::data_out;
			using implicit_plan::matrix_n;
			using implicit_plan::matrix_m;
			using implicit_plan::grid_n;
			using implicit_plan::grid_m;
		
			double alpha; //!< The implicit fraction of the plan (1.0 for purely implicit, 0.0 for purely explicit)
			std::vector <double> new_matrix_vec; //!< An array that contains the diffusion matrix
			double *new_matrix; //!< A pointer to the above vector for convenience
			const double *pos_m; //!< A pointer to the vertical position, for speed
			const double *pos_n; //!< A pointer to the horizontal positions, for speed

		public:
			using implicit_plan::element_flags;
			using implicit_plan::component_flags;
		
			/*!**********************************************************************
			 * \param i_data_in A reference to the input data variable for the plan
			 * \param i_data_out A reference to the output data variable for the plan
			 * \param i_matrix_n The double matrix associated with the horizontal solve
			 * \param i_matrix_m The double matrix associated with the vertical solve
			 * \param i_coeff The diffusion coefficient
			 * \param i_alpha The implicit fraction (1.0 for purely implicit, 0.0 for purely explicit)
			 ************************************************************************/
			vertical (double i_alpha, double *i_matrix_n, double *i_matrix_m, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0) : implicit_plan (i_matrix_n, i_matrix_m, i_data_in, i_data_out, i_coeff, spectral_spectral), alpha (i_alpha) {
				new_matrix_vec.resize (m * m, 0.0);
				new_matrix = &new_matrix_vec [0];

				pos_m = &(grid_m [0]);
				pos_n = &(grid_n [0]);
				setup ();
			}
		
			virtual ~vertical () {}

			void setup () {
				TRACE ("Setting up with coefficient " << coeff);
				if (matrix_m) {
					linalg::scale (m * m, 0.0, new_matrix);
					for (int j = 0; j < m; ++j) {
						linalg::add_scaled (m, coeff, grid_m.get_data (2) + j, new_matrix + j, m, m);
					}
					linalg::add_scaled (m * m, -1.0 * alpha, new_matrix, matrix_m);
				} else {
					WARN ("No matrix");
				}
			}
			
			/*!**********************************************************************
			 * \copydoc implicit_plan::execute
			 ************************************************************************/
			virtual void execute () {
				TRACE ("Operating...");
				// Depending on the direction of the solve, treat this term as either partially or fully explicit
				if (component_flags & z_solve) {
					if (1.0 - alpha != 0.0) {
						linalg::matrix_matrix_multiply (m, ldn, m, 1.0 - alpha, new_matrix, data_in, 1.0, data_out);
					}
				} else {
					linalg::matrix_matrix_multiply (m, ldn, m, 1.0, new_matrix, data_in, 1.0, data_out);
					// #pragma omp parallel for
					// for (int i = 0; i < n; ++i)
					// {
					// 	for (int j = 1; j < m - 1; ++j)
					// 	{
					// 		data_out [i * m + j] += coeff * ((data_in [i * m + j + 1] - data_in [i * m + j]) / (pos_m [j + 1] - pos_m [j]) - (data_in [i * m + j] - data_in [i * m + j - 1]) / (pos_m [j] - pos_m [j - 1])) / (pos_m [j + 1] - pos_m [j - 1]);
					// 	}
					// }

					// #pragma omp parallel for
					// for (int i = 1; i < m - 1; ++i)
					// {
					// 	for (int j = 0; j < m; ++j)
					// 	{
					// 		data_out [i * m + j] += coeff * ((data_in [(i + 1) * m + j] - data_in [i * m + j]) / (pos_n [i + 1] - pos_n [i]) - (data_in [i * m + j] - data_in [(i - 1) * m + j]) / (pos_n [i] - pos_n [i - 1])) / (pos_n [i + 1] - pos_n [i - 1]);
					// 	}
					// }

					// for (int j = 0; j < m; ++j)
					// {
					// 	data_out [j] += coeff * ((data_in [m + j] - data_in [j]) / (pos_n [1] - pos_n [0]) - (data_in [j] - data_in [(n - 1) * m + j]) / (pos_n [1] - pos_n [0])) / 2. / (pos_n [1] - pos_n [0]);
					// 	data_out [(n - 1) * m + j] += coeff * ((data_in [j] - data_in [(n - 1) * m + j]) / (pos_n [n - 1] - pos_n [n - 2]) - (data_in [(n - 1) * m + j] - data_in [(n - 2) * m + j]) / (pos_n [n - 1] - pos_n [n - 2])) / 2. / (pos_n [n - 1] - pos_n [n - 2]);
					// }
				}


				TRACE ("Operation complete.");
			}
		
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public implicit_plan::factory
			{
			private:
				double alpha; //!< The implicit fraction of the plan to be constructed
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The diffusion coefficient of the plan to be constructed
				 * \param i_alpha The implicit fraction of the plan to be constructed
				 ************************************************************************/
				factory (double i_coeff = 1.0, double i_alpha = 1.0) : implicit_plan::factory (i_coeff), alpha (i_alpha) {}

				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plan > (new vertical (alpha, matrices [0], matrices [1], i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plan > ();
				}
			};
		};
	
		/*!**********************************************************************
		 * \brief An implicit plan for vertical diffusion with z-dependence
		 ************************************************************************/
		class background_vertical : public implicit_plan
		{
		protected:
			using implicit_plan::coeff;
			using implicit_plan::n;
			using implicit_plan::ldn;
			using implicit_plan::m;
			using implicit_plan::data_in;
			using implicit_plan::data_out;
			using implicit_plan::matrix_n;
			using implicit_plan::matrix_m;
			using implicit_plan::grid_n;
			using implicit_plan::grid_m;
		
			double alpha; //!< The implicit fraction of the plan (1.0 for purely implicit, 0.0 for purely explicit)
			double *diffusion; //!< A pointer to a vector of diffusion coefficients

			std::vector <double> coeff_dz_vec; //!< A vector of derivatives of the coefficient
			std::vector <double> new_matrix_vec; //!< A vector to contain the diffusion matrix
			double *oodz; //!< One over the boundary centered differential of z
			double *coeff_dz; //!< A pointer to coeff_dz vector
			double *new_matrix; //!< A pointer to the new_matrix vector
		
		public:
			using implicit_plan::element_flags;
			using implicit_plan::component_flags;
		
			/*!**********************************************************************
			 * \param i_data_in A reference to the input data variable for the plan
			 * \param i_data_out A reference to the output data variable for the plan
			 * \param i_matrix_n The double matrix associated with the horizontal solve
			 * \param i_matrix_m The double matrix associated with the vertical solve
			 * \param i_alpha The implicit fraction of the plan (1.0 for purely implicit, 0.0 for purely explicit)
			 * \param i_diffusion A pointer to a vector of diffusion coefficients
			 ************************************************************************/
			background_vertical (double i_alpha, double *i_diffusion, double *i_matrix_n, double *i_matrix_m, grids::variable &i_data_in, grids::variable &i_data_out) : implicit_plan (i_matrix_n, i_matrix_m, i_data_in, i_data_out, 1.0), alpha (i_alpha), diffusion (i_diffusion) {
				oodz = grid_m.get_ood2 ();
				coeff_dz_vec.resize (m);
				coeff_dz = &coeff_dz_vec [0];
				new_matrix_vec.resize (m * m);
				new_matrix = &new_matrix_vec [0];

				setup ();
			}
		
			virtual ~background_vertical () {}

			/**
			 * @copydoc plan::setup
			 */
			void setup () {
				TRACE ("Setting up");
				coeff_dz [0] = coeff * (diffusion [1] - diffusion [0]) * oodz [0];
				for (int i = 1; i < m - 1; ++i) {
					coeff_dz [i] = coeff * (diffusion [i + 1] - diffusion [i - 1]) * oodz [i];
				}
				coeff_dz [m - 1] = coeff * (diffusion [m - 1] - diffusion [m - 2]) * oodz [m - 1];
				
				if (matrix_m) {
					for (int j = 0; j < m; ++j) {
						// DEBUG ("Updating diff " << diffusion [j]);
						linalg::add_scaled (m, coeff * diffusion [j], grid_m.get_data (2) + j, new_matrix + j, m, m);
						linalg::add_scaled (m, coeff_dz [j], grid_m.get_data (1) + j, new_matrix + j, m, m);
					}
					linalg::add_scaled (m * m, -1.0 * alpha, new_matrix, matrix_m);
				} else {
					WARN ("No matrix");
				}
			}
			
			/*!**********************************************************************
			 * \copydoc implicit_plan::execute
			 ************************************************************************/
			void execute () {
				TRACE ("Operating...");
				// Depending on the direction of the solve, treat this term as either partially or fully explicit
				if (component_flags & z_solve) {
					if (1.0 - alpha != 0.0) {
						linalg::matrix_matrix_multiply (m, ldn, m, 1.0 - alpha, new_matrix, data_in, 1.0, data_out, m);
					}
				} else {
					linalg::matrix_matrix_multiply (m, ldn, m, 1.0, new_matrix, data_in, 1.0, data_out, m);
				}

				TRACE ("Operation complete.");
			}
		
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public implicit_plan::factory
			{
			private:
				double alpha; //!< The implicit fraction for the plan to be constructed
				double *diffusion; //!< A pointer to the vector of diffusion coefficients
			
			public:
				/*!**********************************************************************
				 * \param i_alpha The implicit fraction for the plan to be constructed
				 * \param i_diffusion A pointer to the vector of diffusion coefficients
				 ************************************************************************/
				factory (double i_alpha, double *i_diffusion) : alpha (i_alpha), diffusion (i_diffusion) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					return std::shared_ptr <plan > (new background_vertical (alpha, diffusion, matrices [0], matrices [1], i_data_in, i_data_out));
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to add a source term to an equation
		 ************************************************************************/
		class vertical_stress : public real_plan
		{
		private:
			using real_plan::coeff;
			using real_plan::n;
			using real_plan::ldn;
			using real_plan::m;
			using real_plan::dims;
			using real_plan::data_in;
			using real_plan::data_out;
			using real_plan::grid_m;
			using real_plan::grid_n;

			double *density; //!< A pointer to the density data
			double *data_other; //!< A pointer to the other component of the velocity
			double pioL; //!< A shorthand for 2*pi divided by the length of the system
			const double *pos_m; //!< A pointer to the vertical position data
			const double *pos_n; //!< A pointer to the horizontal position data
			double *oodx; //!< An array of one over the differential of x, centered in the cell
			double *oodx2; //!< An array of one over the differential of x, centered at the boundary
			double *oodz; //!< An array of one over the differential of z, centered in the cell
			double *oodz2; //!< An array of one over the differential of z, centered at the boundary
		
		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * @param i_density A reference to the density data
			 * \param i_data_other A reference to the other velocity component
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			vertical_stress (grids::variable &i_density, grids::variable &i_data_other, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0) : 
			real_plan (i_data_in, i_data_out, i_coeff / 3.0), 
			density (i_density.ptr ()),
			data_other (i_data_other.ptr ()) {
				TRACE ("Adding stress...");
				pioL = 2.0 * (std::acos (-1.0) / (grid_n [n - 1] - grid_n [0]));

				pos_n = &grid_n [0];
				pos_m = &grid_m [0];

				oodx = grid_n.get_ood ();
				oodx2 = grid_n.get_ood2 ();
				oodz = grid_m.get_ood ();
				oodz2 = grid_m.get_ood2 ();
			}
		
			virtual ~vertical_stress () {}
		
			/*!**********************************************************************
			 * \copydoc explicit_plan::execute
			 ************************************************************************/
			virtual void execute () {
				TRACE ("Executing source...");
				#pragma omp parallel for
				for (int i = 0; i < n; ++i)
				{
					int p1 = 0, m1 = 0, g = 0;
					p1 = (i + 1) % n;
					m1 = (i - 1 + n) % n;
					for (int j = 1; j < m - 1; ++j)
					{
						g = i * m + j;
						data_out [g] += coeff / 3. * (density [g + 1] + density [g]) / 2. * (data_in [g + 1] - data_in [g]) * oodz [j] * oodz2 [j] / density [g];
						data_out [g] -= coeff / 3. * (density [g] + density [g - 1]) / 2. * (data_in [g] - data_in [g - 1]) * oodz [j - 1] * oodz2 [j] / density [g];

						// data_out [g] -= coeff * 2. / 3. * (density [g + 1] + density [g]) / 2. * (data_other [p1 * m + j] - data_other [g]) * oodx [i] * oodz2 [j] / density [g];
						// data_out [g] += coeff * 2. / 3. * (density [g] + density [g - 1]) / 2. * (data_other [g] - data_other [m1 * m + j]) * oodx [m1] * oodz2 [j] / density [g];

						// data_out [g] += coeff * (density [p1 * m + j] + density [g]) / 2. * (data_other [g + 1] - data_other [g]) * oodz [j] * oodx2 [i] / density [g];
						// data_out [g] -= coeff * (density [g] + density [m1 * m + j]) / 2. * (data_other [g] - data_other [g - 1]) * oodz [j - 1] * oodx2 [i] / density [g];

						data_out [g] -= coeff * 2. / 3. * density [g + 1] * (data_other [p1 * m + j + 1] - data_other [m1 * m + j + 1]) * oodx2 [i] * oodz2 [j] / density [g];
						data_out [g] += coeff * 2. / 3. * density [g - 1] * (data_other [p1 * m + j - 1] - data_other [m1 * m + j - 1]) * oodx2 [i] * oodz2 [j] / density [g];

						data_out [g] += coeff * density [p1 * m + j] * (data_other [p1 * m + j + 1] - data_other [p1 * m + j - 1]) * oodz2 [j] * oodx2 [i] / density [g];
						data_out [g] -= coeff * density [m1 * m + j] * (data_other [m1 * m + j + 1] - data_other [m1 * m + j - 1]) * oodz2 [j] * oodx2 [i] / density [g];
					}
				}
			}
		
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public explicit_plan::factory
			{
			private:
				grids::variable &density; //!< The data density to be used when constructing the plan
				grids::variable &data_other; //!< The other velocity component to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_density A reference to the density of the system
				 * \param i_data_other The other velocity component to be used when constructing the plan
				 ************************************************************************/
				factory (grids::variable &i_density, grids::variable &i_data_other, double i_coeff = 1.0) : 
				real_plan::factory (i_coeff), 
				density (i_density),
				data_other (i_data_other) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan > (new vertical_stress (density, data_other, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan > ();
				}
			};
		};
	} /* diffusion */
} /* plans */

#endif /* end of include guard: DIFFUSION_TWO_D_HPP_YHECX9VS */
