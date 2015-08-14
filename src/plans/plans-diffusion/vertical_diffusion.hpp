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
		template <class datatype>
		class vertical : public implicit_plan <datatype>
		{
		protected:
			using implicit_plan <datatype>::coeff;
			using implicit_plan <datatype>::n;
			using implicit_plan <datatype>::ldn;
			using implicit_plan <datatype>::m;
			using implicit_plan <datatype>::dims;
			using implicit_plan <datatype>::data_in;
			using implicit_plan <datatype>::data_out;
			using implicit_plan <datatype>::matrix_n;
			using implicit_plan <datatype>::matrix_m;
			using implicit_plan <datatype>::grid_n;
			using implicit_plan <datatype>::grid_m;
		
			datatype alpha; //!< The implicit fraction of the plan (1.0 for purely implicit, 0.0 for purely explicit)
			std::vector <datatype> new_matrix_vec;
			datatype *new_matrix;
			const datatype *pos_m;
			const datatype *pos_n;

		public:
			using implicit_plan <datatype>::element_flags;
			using implicit_plan <datatype>::component_flags;
		
			/*!**********************************************************************
			 * \copydoc implicit_plan::implicit_plan
			 * 
			 * \param i_coeff The datatype diffusion coefficient
			 * \param i_alpha The implicit fraction of the plan (1.0 for purely implicit, 0.0 for purely explicit)
			 ************************************************************************/
			vertical (datatype i_alpha, datatype *i_matrix_n, datatype *i_matrix_m, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0) : implicit_plan <datatype> (i_matrix_n, i_matrix_m, i_data_in, i_data_out, i_coeff, spectral_spectral), alpha (i_alpha) {
				new_matrix_vec.resize (m * m, 0.0);
				new_matrix = &new_matrix_vec [0];

				pos_m = &(grid_m [0]);
				pos_n = &(grid_n [0]);
				setup ();
			}
		
			virtual ~vertical () {}

			void setup () {
				TRACE ("Setting up with coefficient " << coeff);
				INFO ("COEFFICIENT VERTICAL IS " << coeff);
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
			 * \copydoc implicit_plan::factory
			 ************************************************************************/
			class factory : public implicit_plan <datatype>::factory
			{
			private:
				datatype alpha; //!< The implicit fraction of the plan to be constructed
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The diffusion coefficient of the plan to be constructed
				 * \param i_alpha The implicit fraction of the plan to be constructed
				 ************************************************************************/
				factory (datatype i_coeff = 1.0, datatype i_alpha = 1.0) : implicit_plan <datatype>::factory (i_coeff), alpha (i_alpha) {}

				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc implicit_plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plan <datatype> > (new vertical <datatype> (alpha, matrices [0], matrices [1], i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plan <datatype> > ();
				}
			};
		};
	
		/*!**********************************************************************
		 * \brief An implicit plan for vertical diffusion with z-dependence
		 ************************************************************************/
		template <class datatype>
		class background_vertical : public implicit_plan <datatype>
		{
		protected:
			using implicit_plan <datatype>::coeff;
			using implicit_plan <datatype>::n;
			using implicit_plan <datatype>::ldn;
			using implicit_plan <datatype>::m;
			using implicit_plan <datatype>::data_in;
			using implicit_plan <datatype>::data_out;
			using implicit_plan <datatype>::matrix_n;
			using implicit_plan <datatype>::matrix_m;
			using implicit_plan <datatype>::grid_n;
			using implicit_plan <datatype>::grid_m;
		
			datatype alpha; //!< The implicit fraction of the plan (1.0 for purely implicit, 0.0 for purely explicit)
			datatype *diffusion; //!< A pointer to a vector of diffusion coefficients
			bool explicit_calculate;

			std::vector <datatype> oodz_vec;
			std::vector <datatype> coeff_dz_vec;
			std::vector <datatype> new_matrix_vec;
			datatype *oodz, *coeff_dz, *new_matrix;
		
		public:
			using implicit_plan <datatype>::element_flags;
			using implicit_plan <datatype>::component_flags;
		
			/*!**********************************************************************
			 * \copydoc implicit_plan::implicit_plan
			 * 
			 * \param i_alpha The implicit fraction of the plan (1.0 for purely implicit, 0.0 for purely explicit)
			 * \param i_diffusion A pointer to a vector of diffusion coefficients
			 ************************************************************************/
			background_vertical (datatype i_alpha, datatype *i_diffusion, bool i_explicit_calculate, datatype *i_matrix_n, datatype *i_matrix_m, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) : implicit_plan <datatype> (i_matrix_n, i_matrix_m, i_data_in, i_data_out, 1.0), alpha (i_alpha), diffusion (i_diffusion), explicit_calculate (i_explicit_calculate) {
				oodz_vec.resize (m);
				oodz = &oodz_vec [0];
				coeff_dz_vec.resize (m);
				coeff_dz = &coeff_dz_vec [0];
				new_matrix_vec.resize (m * m);
				new_matrix = &new_matrix_vec [0];

				oodz [0] = 1.0 / (grid_m [1] - grid_m [0]);
				for (int i = 1; i < m-1; ++i) {
					oodz [i] = 1.0 / (grid_m [i + 1] - grid_m [i - 1]);
				}
				oodz [m - 1] = 1.0 / (grid_m [m - 1] - grid_m [m - 2]);
				
				setup ();
			}
		
			virtual ~background_vertical () {}

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
			 * \copydoc implicit_plan::factory
			 ************************************************************************/
			class factory : public implicit_plan <datatype>::factory
			{
			private:
				datatype alpha; //!< The implicit fraction for the plan to be constructed
				datatype *diffusion; //!< A pointer to the vector of diffusion coefficients
				bool explicit_calculate; //!< A pointer to the vector of diffusion coefficients
			
			public:
				/*!**********************************************************************
				 * \param i_alpha The implicit fraction for the plan to be constructed
				 * \param i_diffusion A pointer to the vector of diffusion coefficients
				 ************************************************************************/
				factory (datatype i_alpha, datatype *i_diffusion, bool i_explicit_calculate = true) : alpha (i_alpha), diffusion (i_diffusion), explicit_calculate (i_explicit_calculate) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc implicit_plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					return std::shared_ptr <plan <datatype> > (new background_vertical <datatype> (alpha, diffusion, explicit_calculate, matrices [0], matrices [1], i_data_in, i_data_out));
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to add a source term to an equation
		 ************************************************************************/
		template <class datatype>
		class vertical_stress : public real_plan <datatype>
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

			datatype *density;
			datatype *data_other;
			datatype pioL;
			const datatype *pos_m, *pos_n;
			datatype *oodx, *oodx2, *oodz, *oodz2;
		
		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * \param i_data_source The data pointer for the source data
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			vertical_stress (grids::variable <datatype> &i_density, grids::variable <datatype> &i_data_other, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0) : 
			real_plan <datatype> (i_data_in, i_data_out, i_coeff / 3.0), 
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
			 * \copydoc explicit_plan::factory
			 ************************************************************************/
			class factory : public explicit_plan <datatype>::factory
			{
			private:
				grids::variable <datatype> &density; //!< The data source to be used when constructing the plan
				grids::variable <datatype> &data_other; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 ************************************************************************/
				factory (grids::variable <datatype> &i_density, grids::variable <datatype> &i_data_other, datatype i_coeff = 1.0) : 
				explicit_plan <datatype>::factory (i_coeff), 
				density (i_density),
				data_other (i_data_other) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc explicit_plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan <datatype> > (new vertical_stress <datatype> (density, data_other, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan <datatype> > ();
				}
			};
		};
	} /* diffusion */
} /* plans */

#endif /* end of include guard: DIFFUSION_TWO_D_HPP_YHECX9VS */
