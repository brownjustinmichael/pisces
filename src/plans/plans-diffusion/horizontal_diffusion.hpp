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
#include "../explicit_plan.hpp"
#include "../real_plan.hpp"
#include "io/parameters.hpp"

namespace plans
{
	namespace diffusion
	{
		/*!**********************************************************************
		 * \brief An implicit plan that calculates horizontal diffusion for Fourier modes
		 ************************************************************************/
		class horizontal : public implicit_plan
		{
		private:
			using implicit_plan::coeff;
			using implicit_plan::n;
			using implicit_plan::ldn;
			using implicit_plan::m;
			using implicit_plan::dims;
			using implicit_plan::grid_n;
			using implicit_plan::data_in;
			using implicit_plan::data_out;
			using implicit_plan::matrix_n;
		
			double alpha; //!< The implicit fraction (1.0 for purely implicit, 0.0 for purely explicit)
			double pioL2; //!< The value pi / L ^ 2 to save computational time
		
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
			horizontal (double i_alpha, double *i_matrix_n, double *i_matrix_m, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0) : implicit_plan (i_matrix_n, i_matrix_m, i_data_in, i_data_out, i_coeff, real_spectral), alpha (i_alpha) {
				TRACE ("Instantiating...");
				pioL2 = 4.0 * (std::acos (-1.0) * std::acos (-1.0) / (grid_n [n - 1] - grid_n [0]) / (grid_n [n - 1] - grid_n [0]));
				setup ();
			}
			
			virtual ~horizontal () {}

			void setup () {
				TRACE ("Setting up");
				if (matrix_n) {
					// For Fourier modes, the matrix is diagonal and not particularly complicated
					// We set up m of these matrices in case there is some z-dependence added in later
					for (int j = 0; j < m * dims; ++j) {
						matrix_n [j] = matrix_n [m + j] = 0.0;
						for (int i = 2; i < ldn; ++i) {
							matrix_n [i * m + j] = coeff * alpha * pioL2 * (double) ((i / 2) * (i / 2));
						}
					}
				} else {
					WARN ("No matrix");
				}
			}
			
			/*!**********************************************************************
			 * \copydoc implicit_plan::execute
			 ************************************************************************/
			void execute () {	
				TRACE ("Operating..." << element_flags);
				// Depending on the direction of the solve, we'll either treat the term as partially implicit or fully explicit
				if (component_flags & x_solve) {
					if (1.0 - alpha != 0.0) {
						// #pragma omp parallel for
						for (int i = 2; i < ldn; ++i) {
							linalg::add_scaled (m * dims, -coeff * (1.0 - alpha) * pioL2 * (i / 2) * (i / 2), data_in + i * m * dims, data_out + i * m * dims);
						}
					}
				} else {
					// #pragma omp parallel for
					for (int i = 2; i < ldn; ++i) {
						linalg::add_scaled (m * dims, -coeff * pioL2 * (i / 2) * (i / 2), data_in + i * m * dims, data_out + i * m * dims);
					}
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
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The diffusion coefficient for the plan to be constructed
				 * \param i_alpha The implicit fraction for the plan to be constructed
				 ************************************************************************/
				factory (double i_coeff = 1.0, double i_alpha = 1.0) : implicit_plan::factory (i_coeff), alpha (i_alpha) {
				}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plan > (new horizontal (alpha, matrices [0], matrices [1], i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plan > ();
				}
			};
		};
	
		/*!**********************************************************************
		 * \brief A horizontal diffusion plan that varies with z
		 ************************************************************************/
		class background_horizontal : public implicit_plan
		{
		private:
			using implicit_plan::coeff;
			using implicit_plan::n;
			using implicit_plan::ldn;
			using implicit_plan::m;
			using implicit_plan::dims;
			using implicit_plan::grid_n;
			using implicit_plan::data_in;
			using implicit_plan::data_out;
			using implicit_plan::matrix_n;
		
			double alpha; //!< The implicit fraction (1.0 for purely implicit, 0.0 for purely explicit)
			double *diffusion; //!< A pointer to the double vector of diffusion coefficients
			double pioL2; //!< The value pi / L ^ 2 to save computational time

		public:
			using implicit_plan::element_flags;
			using implicit_plan::component_flags;
		
			/*!**********************************************************************
			 * \param i_data_in A reference to the input data variable for the plan
			 * \param i_data_out A reference to the output data variable for the plan
			 * \param i_matrix_n The double matrix associated with the horizontal solve
			 * \param i_matrix_m The double matrix associated with the vertical solve
			 * \param i_alpha The implicit fraction (1.0 for purely implicit, 0.0 for purely explicit)
			 * \param i_diffusion A pointer to the double vector of diffusion coefficients
			 ************************************************************************/
			background_horizontal (double i_alpha, double *i_diffusion, double *i_matrix_n, double *i_matrix_m, grids::variable &i_data_in, grids::variable &i_data_out) : implicit_plan (i_matrix_n, i_matrix_m, i_data_in, i_data_out, 1.0, real_spectral), alpha (i_alpha), diffusion (i_diffusion) {
				TRACE ("Instantiating...");
				pioL2 = 4.0 * (std::acos (-1.0) * std::acos (-1.0) / (grid_n [n - 1] - grid_n [0]) / (grid_n [n - 1] - grid_n [0]));
				setup ();
			}
			
			virtual ~background_horizontal () {}
			
			void setup () {
				// For Fourier modes, the matrix is diagonal and not particularly complicated
				TRACE ("Setting up");
				if (matrix_n) {
					for (int j = 0; j < m; ++j) {
						for (int k = 0; k < dims; ++k)
						{
							matrix_n [j] = matrix_n [(m + j) * dims + k] = 0.0;
							for (int i = 2; i < ldn; ++i) {
								matrix_n [(i * m + j) * dims + k] = coeff * diffusion [j] * alpha * pioL2 * (double) ((i / 2) * (i / 2));
							}						
						}
					}
				} else {
					WARN ("No matrix");
				}
			}

			/*!**********************************************************************
			 * \copydoc implicit_plan::execute
			 ************************************************************************/
			void execute () {	
				TRACE ("Operating..." << element_flags);
				// Depending on the direction of the solve, we'll either treat the term as partially implicit or fully explicit
				if (component_flags & x_solve) {
					if (1.0 - alpha != 0.0) {
						// #pragma omp parallel for
						for (int j = 0; j < m; ++j) {
							for (int k = 0; k < dims; ++k)
							{
								for (int i = 2; i < ldn; ++i) {
									data_out [(i * m + j) * dims + k] -= coeff * diffusion [j] * (1.0 - alpha) * pioL2 * (i / 2) * (i / 2) * data_in [(i * m + j) * dims + k];
								}
							}
						}
					}
				} else {
					// #pragma omp parallel for
					for (int j = 0; j < m; ++j) {
						for (int i = 2; i < ldn; ++i) {
							for (int k = 0; k < dims; ++k)
							{
								data_out [(i * m + j) * dims + k] -= coeff * diffusion [j] * pioL2 * (i / 2) * (i / 2) * data_in [(i * m + j) * dims + k];
							}
						}
					}
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
				double *diffusion; //!< The pointer to the diffusion vector for the plan to be constructed
			
			public:
				/*!**********************************************************************
				 * \param i_alpha The implicit fraction for the plan to be constructed
				 * \param i_diffusion The pointer to the diffusion vector for the plan to be constructed
				 ************************************************************************/
				factory (double i_alpha, double *i_diffusion) : alpha (i_alpha), diffusion (i_diffusion) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					return std::shared_ptr <plan > (new background_horizontal (alpha, diffusion, matrices [0], matrices [1], i_data_in, i_data_out));
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to add a source term to an equation
		 ************************************************************************/
		class horizontal_stress : public real_plan
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
			double *data_other; //!< A pointer to the other velocity component
			double pioL; //!< The value 2*pi / the width of the system, for speed and convenience
			const double *pos_m; //!< A pointer to the vertical position, for speed
			const double *pos_n; //!< A pointer to the horizontal positions, for speed
			double *oodx; //!< An array of one over the differential of x, centered in the cell
			double *oodx2; //!< An array of one over the differential of x, centered at the boundary
			double *oodz; //!< An array of one over the differential of z, centered in the cell
			double *oodz2; //!< An array of one over the differential of z, centered at the boundary
		
		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * \param i_data_other The data pointer for the other velocity component
			 * @param i_density A reference to the density variable
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			horizontal_stress (grids::variable &i_density, grids::variable &i_data_other, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0) : 
			real_plan (i_data_in, i_data_out, i_coeff), 
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
		
			virtual ~horizontal_stress () {}

			/*!**********************************************************************
			 * \copydoc explicit_plan::execute
			 ************************************************************************/
			virtual void execute () {
				TRACE ("Executing source...");
				DEBUG ("The pointer is now " << data_other);
				#pragma omp parallel for
				for (int i = 0; i < n; ++i)
				{
					int p1 = 0, m1 = 0, g = 0;
					p1 = (i + 1) % n;
					m1 = (i - 1 + n) % n;
					for (int j = 1; j < m - 1; ++j)
					{
						g = i * m + j;
						data_out [g] += coeff / 3. * (density [p1 * m + j] + density [g]) / 2. * (data_in [p1 * m + j] - data_in [g]) * oodx [i] * oodx2 [i] / density [g];
						data_out [g] -= coeff / 3. * (density [g] + density [m1 * m + j]) / 2. * (data_in [g] - data_in [m1 * m + j]) * oodx [m1] * oodx2 [i] / density [g];

						data_out [g] -= coeff * 2. / 3. * density [p1 * m + j] * (data_other [p1 * m + j + 1] - data_other [p1 * m + j - 1]) * oodz2 [j] * oodx2 [i] / density [g];
						data_out [g] += coeff * 2. / 3. * density [m1 * m + j] * (data_other [m1 * m + j + 1] - data_other [m1 * m + j - 1]) * oodz2 [j] * oodx2 [i] / density [g];

						data_out [g] += coeff * density [g + 1] * (data_other [p1 * m + j + 1] - data_other [m1 * m + j + 1]) * oodx2 [i] * oodz2 [j] / density [g];
						data_out [g] -= coeff * density [g - 1] * (data_other [p1 * m + j - 1] - data_other [m1 * m + j - 1]) * oodx2 [i] * oodz2 [j] / density [g];
					}
				}
			}
		
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public real_plan::factory
			{
			private:
				grids::variable &density; //!< The data source to be used when constructing the plan
				grids::variable &data_other; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_other The data of the other velocity component to be used when constructing the plan
				 * @param i_density A reference to the density variable
				 ************************************************************************/
				factory (grids::variable &i_density, grids::variable &i_data_other, double i_coeff = 1.0) : 
				explicit_plan::factory (i_coeff), 
				density (i_density),
				data_other (i_data_other) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan > (new horizontal_stress (density, data_other, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan > ();
				}
			};
		};
	} /* diffusion */
} /* plans */

#endif /* end of include guard: HORIZONTAL_DIFFUSION_HPP_E1A9847D */
