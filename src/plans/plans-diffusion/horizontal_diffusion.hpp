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
			 * \copydoc implicit_plan::implicit_plan
			 * 
			 * \param i_coeff The diffusion coefficient
			 * \param i_alpha The implicit fraction (1.0 for purely implicit, 0.0 for purely explicit)
			 ************************************************************************/
			horizontal (datatype i_alpha, datatype *i_matrix_n, datatype *i_matrix_m, grids::variable <datatype> &i_data_in, datatype *i_data_out = NULL, datatype i_coeff = 1.0, int *i_element_flags = NULL, int *i_component_flags = NULL) : implicit_plan <datatype> (i_matrix_n, i_matrix_m, i_data_in, i_data_out, i_coeff, i_element_flags, i_component_flags), alpha (i_alpha) {
				TRACE ("Instantiating...");
				pioL2 = 4.0 * (std::acos (-1.0) * std::acos (-1.0) / (grid_n [n - 1] - grid_n [0]) / (grid_n [n - 1] - grid_n [0]));
				setup ();
			}
			
			virtual ~horizontal () {}
			
			virtual int type () {
				return plan <datatype>::mid;
			}

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
				if (*component_flags & x_solve) {
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
			 * \copydoc implicit_plan::factory
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
				factory (datatype i_coeff = 1.0, datatype i_alpha = 1.0) : implicit_plan <datatype>::factory (i_coeff), alpha (i_alpha) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc implicit::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
					if (coeff) {
						return std::shared_ptr <plan <datatype> > (new horizontal <datatype> (alpha, matrices [0], matrices [1], i_data_in, i_data_out, 1.0, i_element_flags, i_component_flags));
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
			 * \copydoc implicit_plan::implicit_plan
			 * 
			 * \param i_alpha The implicit fraction (1.0 for purely implicit, 0.0 for purely explicit)
			 * \param i_diffusion A pointer to the datatype vector of diffusion coefficients
			 ************************************************************************/
			background_horizontal (datatype i_alpha, datatype *i_diffusion, datatype *i_matrix_n, datatype *i_matrix_m, grids::variable <datatype> &i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) : implicit_plan <datatype> (i_matrix_n, i_matrix_m, i_data_in, i_data_out, 1.0, i_element_flags, i_component_flags), alpha (i_alpha), diffusion (i_diffusion) {
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

			virtual int type () {
				return plan <datatype>::mid;
			}
			
			/*!**********************************************************************
			 * \copydoc implicit_plan::execute
			 ************************************************************************/
			void execute () {	
				TRACE ("Operating..." << element_flags);
				// Depending on the direction of the solve, we'll either treat the term as partially implicit or fully explicit
				if (*component_flags & x_solve) {
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
			 * \copydoc implicit_plan::factory
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
				 * \copydoc implicit::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
					return std::shared_ptr <plan <datatype> > (new background_horizontal <datatype> (alpha, diffusion, matrices [0], matrices [1], i_data_in, i_data_out, i_element_flags, i_component_flags));
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to add a source term to an equation
		 ************************************************************************/
		template <class datatype>
		class horizontal_stress : public explicit_plan <datatype>
		{
		private:
			using explicit_plan <datatype>::coeff;
			using explicit_plan <datatype>::n;
			using explicit_plan <datatype>::ldn;
			using explicit_plan <datatype>::m;
			using explicit_plan <datatype>::dims;
			using explicit_plan <datatype>::data_in;
			using explicit_plan <datatype>::data_out;
			using explicit_plan <datatype>::grid_m;
			using explicit_plan <datatype>::grid_n;

			datatype *data_other;
			datatype pioL;
			const datatype *pos_m;
			std::vector <datatype> diff;
			std::vector <datatype> diff2;
		
		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * \param i_data_source The data pointer for the source data
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			horizontal_stress (grids::variable <datatype> &i_data_other, grids::variable <datatype> &i_data_in, datatype *i_data_out, datatype i_coeff = 1.0, int *i_element_flags = NULL, int *i_component_flags = NULL) : 
			explicit_plan <datatype> (i_data_in, i_data_out, i_element_flags, i_component_flags, i_coeff), 
			data_other (i_data_other.ptr ()) {
				TRACE ("Adding stress...");
				pioL = 2.0 * (std::acos (-1.0) / (grid_n [n - 1] - grid_n [0]));

				pos_m = &grid_m [0];

				diff.resize (m);
				for (int j = 0; j < m - 1; ++j)
				{
					diff [j] = pos_m [j + 1] - pos_m [j];
				}
				diff [m - 1] = pos_m [m - 1] - pos_m [m - 2];

				diff2.resize (m);
				for (int j = 1; j < m - 1; ++j)
				{
					diff2 [j] = pos_m [j + 1] - pos_m [j - 1];
				}
				diff2 [0] = pos_m [1] - pos_m [0];
				diff2 [m - 1] = pos_m [m - 1] - pos_m [m - 2];
			}
		
			virtual ~horizontal_stress () {}

			virtual int type () {
				return plan <datatype>::mid;
			}
		
			/*!**********************************************************************
			 * \copydoc explicit_plan::execute
			 ************************************************************************/
			virtual void execute () {
				TRACE ("Executing source...");
				for (int i = 0; i < ldn; i += 2)
				{
					for (int j = 1; j < m - 1; ++j)
					{
						data_out [i * m + j] += coeff / 3. * pioL * pioL * (i / 2) * (i / 2) * data_in [i * m + j];
						data_out [i * m + j] -= coeff / 3. * pioL * (data_other [(i + 1) * m + j + 1] - data_other [(i + 1) * m + j - 1]) / diff2 [j];

						data_out [(i + 1) * m + j] += coeff / 3. * pioL * pioL * (i / 2) * (i / 2) * data_in [(i + 1) * m + j];
						data_out [(i + 1) * m + j] += coeff / 3. * pioL * (data_other [i * m + j + 1] - data_other [i * m + j - 1]) / diff2 [j];
					}
					// data_out [i * m + j] += coeff / 3. * pioL * pioL * (i / 2) * (i / 2) * data_in [i * m];
					// data_out [i * m + j] -= coeff / 3. * pioL * (data_other [(i + 1) * m + j + 1] - data_other [(i + 1) * m + j - 1]) / diff2 [0];
					
					// data_out [(i + 1) * m + j] += coeff / 3. * pioL * pioL * (i / 2) * (i / 2) * data_in [(i + 1) * m + j + 1];
					// data_out [(i + 1) * m + j] += coeff / 3. * pioL * (data_other [i * m + j + 1] - data_other [i * m + j - 1]) / diff2 [j];
				}
			}
		
			/*!**********************************************************************
			 * \copydoc explicit_plan::factory
			 ************************************************************************/
			class factory : public explicit_plan <datatype>::factory
			{
			private:
				grids::variable <datatype> &data_other; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 ************************************************************************/
				factory (grids::variable <datatype> &i_data_other, datatype i_coeff = 1.0) : 
				explicit_plan <datatype>::factory (i_coeff), 
				data_other (i_data_other) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc explicit_plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
					if (coeff) {
						return std::shared_ptr <plans::plan <datatype> > (new horizontal_stress <datatype> (data_other, i_data_in, i_data_out, 1.0, i_element_flags, i_component_flags));
					}
					return std::shared_ptr <plans::plan <datatype> > ();
				}
			};
		};
	} /* diffusion */
} /* plans */

#endif /* end of include guard: HORIZONTAL_DIFFUSION_HPP_E1A9847D */
