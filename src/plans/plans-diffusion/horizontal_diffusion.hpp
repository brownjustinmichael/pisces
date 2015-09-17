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
	namespace diffusion
	{
		/*!**********************************************************************
		 * \brief An implicit plan that calculates horizontal diffusion for Fourier modes
		 ************************************************************************/
		template <class datatype>
		class horizontal : public implicit_plan <datatype>
		{
		private:
			using implicit_plan <datatype>::coeff;
			using implicit_plan <datatype>::n;
			using implicit_plan <datatype>::ldn;
			using implicit_plan <datatype>::m;
			using implicit_plan <datatype>::dims;
			using implicit_plan <datatype>::grid_n;
			using implicit_plan <datatype>::data_in;
			using implicit_plan <datatype>::data_out;
			using implicit_plan <datatype>::matrix_n;
		
			datatype alpha; //!< The implicit fraction (1.0 for purely implicit, 0.0 for purely explicit)
			datatype pioL2; //!< The value pi / L ^ 2 to save computational time
		
		public:
			using implicit_plan <datatype>::element_flags;
			using implicit_plan <datatype>::component_flags;
		
			/*!**********************************************************************
			 * \param i_data_in A reference to the input data variable for the plan
			 * \param i_data_out A reference to the output data variable for the plan
			 * \param i_matrix_n The datatype matrix associated with the horizontal solve
			 * \param i_matrix_m The datatype matrix associated with the vertical solve
			 * \param i_coeff The diffusion coefficient
			 * \param i_alpha The implicit fraction (1.0 for purely implicit, 0.0 for purely explicit)
			 ************************************************************************/
			horizontal (datatype i_alpha, datatype *i_matrix_n, datatype *i_matrix_m, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0) : implicit_plan <datatype> (i_matrix_n, i_matrix_m, i_data_in, i_data_out, i_coeff, real_spectral), alpha (i_alpha) {
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
							matrix_n [i * m + j] = coeff * alpha * pioL2 * (datatype) ((i / 2) * (i / 2));
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
			class factory : public implicit_plan <datatype>::factory
			{
			private:
				datatype alpha; //!< The implicit fraction for the plan to be constructed
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The diffusion coefficient for the plan to be constructed
				 * \param i_alpha The implicit fraction for the plan to be constructed
				 ************************************************************************/
				factory (datatype i_coeff = 1.0, datatype i_alpha = 1.0) : implicit_plan <datatype>::factory (i_coeff), alpha (i_alpha) {
				}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plan <datatype> > (new horizontal <datatype> (alpha, matrices [0], matrices [1], i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plan <datatype> > ();
				}
			};
		};
	
		/*!**********************************************************************
		 * \brief A horizontal diffusion plan that varies with z
		 ************************************************************************/
		template <class datatype>
		class background_horizontal : public implicit_plan <datatype>
		{
		private:
			using implicit_plan <datatype>::coeff;
			using implicit_plan <datatype>::n;
			using implicit_plan <datatype>::ldn;
			using implicit_plan <datatype>::m;
			using implicit_plan <datatype>::dims;
			using implicit_plan <datatype>::grid_n;
			using implicit_plan <datatype>::data_in;
			using implicit_plan <datatype>::data_out;
			using implicit_plan <datatype>::matrix_n;
		
			datatype alpha; //!< The implicit fraction (1.0 for purely implicit, 0.0 for purely explicit)
			datatype *diffusion; //!< A pointer to the datatype vector of diffusion coefficients
			datatype pioL2; //!< The value pi / L ^ 2 to save computational time

		public:
			using implicit_plan <datatype>::element_flags;
			using implicit_plan <datatype>::component_flags;
		
			/*!**********************************************************************
			 * \param i_data_in A reference to the input data variable for the plan
			 * \param i_data_out A reference to the output data variable for the plan
			 * \param i_matrix_n The datatype matrix associated with the horizontal solve
			 * \param i_matrix_m The datatype matrix associated with the vertical solve
			 * \param i_alpha The implicit fraction (1.0 for purely implicit, 0.0 for purely explicit)
			 * \param i_diffusion A pointer to the datatype vector of diffusion coefficients
			 ************************************************************************/
			background_horizontal (datatype i_alpha, datatype *i_diffusion, datatype *i_matrix_n, datatype *i_matrix_m, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) : implicit_plan <datatype> (i_matrix_n, i_matrix_m, i_data_in, i_data_out, 1.0, real_spectral), alpha (i_alpha), diffusion (i_diffusion) {
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
								matrix_n [(i * m + j) * dims + k] = coeff * diffusion [j] * alpha * pioL2 * (datatype) ((i / 2) * (i / 2));
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
			class factory : public implicit_plan <datatype>::factory
			{
			private:
				datatype alpha; //!< The implicit fraction for the plan to be constructed
				datatype *diffusion; //!< The pointer to the diffusion vector for the plan to be constructed
			
			public:
				/*!**********************************************************************
				 * \param i_alpha The implicit fraction for the plan to be constructed
				 * \param i_diffusion The pointer to the diffusion vector for the plan to be constructed
				 ************************************************************************/
				factory (datatype i_alpha, datatype *i_diffusion) : alpha (i_alpha), diffusion (i_diffusion) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					return std::shared_ptr <plan <datatype> > (new background_horizontal <datatype> (alpha, diffusion, matrices [0], matrices [1], i_data_in, i_data_out));
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to add a source term to an equation
		 ************************************************************************/
		template <class datatype>
		class horizontal_stress : public real_plan <datatype>
		{
		private:
			using real_plan <datatype>::coeff;
			using real_plan <datatype>::n;
			using real_plan <datatype>::ldn;
			using real_plan <datatype>::m;
			using real_plan <datatype>::dims;
			using real_plan <datatype>::data_in;
			using real_plan <datatype>::data_out;
			using real_plan <datatype>::grid_m;
			using real_plan <datatype>::grid_n;

			datatype *density; //!< A pointer to the density data
			datatype *data_other; //!< A pointer to the other velocity component
			datatype pioL; //!< The value 2*pi / the width of the system, for speed and convenience
			const datatype *pos_m; //!< A pointer to the vertical position, for speed
			const datatype *pos_n; //!< A pointer to the horizontal positions, for speed
			datatype *oodx; //!< An array of one over the differential of x, centered in the cell
			datatype *oodx2; //!< An array of one over the differential of x, centered at the boundary
			datatype *oodz; //!< An array of one over the differential of z, centered in the cell
			datatype *oodz2; //!< An array of one over the differential of z, centered at the boundary
		
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
			horizontal_stress (grids::variable <datatype> &i_density, grids::variable <datatype> &i_data_other, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0) : 
			real_plan <datatype> (i_data_in, i_data_out, i_coeff), 
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
			class factory : public explicit_plan <datatype>::factory
			{
			private:
				grids::variable <datatype> &density; //!< The data source to be used when constructing the plan
				grids::variable <datatype> &data_other; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_other The data of the other velocity component to be used when constructing the plan
				 * @param i_density A reference to the density variable
				 ************************************************************************/
				factory (grids::variable <datatype> &i_density, grids::variable <datatype> &i_data_other, datatype i_coeff = 1.0) : 
				explicit_plan <datatype>::factory (i_coeff), 
				density (i_density),
				data_other (i_data_other) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan <datatype> > (new horizontal_stress <datatype> (density, data_other, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan <datatype> > ();
				}
			};
		};
	} /* diffusion */
} /* plans */

#endif /* end of include guard: HORIZONTAL_DIFFUSION_HPP_E1A9847D */
