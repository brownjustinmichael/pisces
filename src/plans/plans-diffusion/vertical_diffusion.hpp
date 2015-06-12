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
		private:
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
		
		public:
			using implicit_plan <datatype>::element_flags;
			using implicit_plan <datatype>::component_flags;
		
			/*!**********************************************************************
			 * \copydoc implicit_plan::implicit_plan
			 * 
			 * \param i_coeff The datatype diffusion coefficient
			 * \param i_alpha The implicit fraction of the plan (1.0 for purely implicit, 0.0 for purely explicit)
			 ************************************************************************/
			vertical (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype i_alpha, datatype *i_matrix_n, datatype *i_matrix_m, datatype *i_data_in, datatype *i_data_out = NULL, datatype i_coeff = 1.0, int *i_element_flags = NULL, int *i_component_flags = NULL) : implicit_plan <datatype> (i_grid_n, i_grid_m, i_matrix_n, i_matrix_m, i_data_in, i_data_out, i_coeff, i_element_flags, i_component_flags), alpha (i_alpha) {
				setup ();
			}
		
			virtual ~vertical () {}

			virtual int type () {
				return plan <datatype>::pre;
			}
			
			void setup () {
				TRACE ("Setting up");
				if (matrix_m) {
					for (int j = 0; j < m; ++j) {
						linalg::add_scaled (m, -coeff * alpha, grid_m.get_data (2) + j, matrix_m + j, m, m);
					}
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
				if (*component_flags & z_solve) {
					if (1.0 - alpha != 0.0) {
						linalg::matrix_matrix_multiply (m, ldn, m, coeff * (1.0 - alpha), grid_m.get_data (2), data_in, 1.0, data_out, m);
					}
				} else {
					linalg::matrix_matrix_multiply (m, ldn, m, coeff, grid_m.get_data (2), data_in, 1.0, data_out, m);
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
				datatype *data_in;
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The diffusion coefficient of the plan to be constructed
				 * \param i_alpha The implicit fraction of the plan to be constructed
				 ************************************************************************/
				factory (datatype *i_data_in, datatype i_coeff = 1.0, datatype i_alpha = 1.0) : implicit_plan <datatype>::factory (i_coeff), alpha (i_alpha), data_in (i_data_in) {}

				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc implicit_plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plan <datatype> > _instance (grids::grid <datatype> **grids, datatype **matrices, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
					DEBUG ("Newer data, " << data_in);

					if (coeff) {
						return std::shared_ptr <plan <datatype> > (new vertical <datatype> (*grids [0], *grids [1], alpha, matrices [0], matrices [1], data_in, i_data_out, 1.0, i_element_flags, i_component_flags));
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
		private:
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
			background_vertical (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype i_alpha, datatype *i_diffusion, bool i_explicit_calculate, datatype *i_matrix_n, datatype *i_matrix_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) : implicit_plan <datatype> (i_grid_n, i_grid_m, i_matrix_n, i_matrix_m, i_data_in, i_data_out, 1.0, i_element_flags, i_component_flags), alpha (i_alpha), diffusion (i_diffusion), explicit_calculate (i_explicit_calculate) {
				oodz_vec.resize (m);
				oodz = &oodz_vec [0];
				coeff_dz_vec.resize (m);
				coeff_dz = &coeff_dz_vec [0];
				new_matrix_vec.resize (m * m);
				new_matrix = &new_matrix_vec [0];

				oodz [0] = 1.0 / (i_grid_m [1] - i_grid_m [0]);
				for (int i = 1; i < m-1; ++i) {
					oodz [i] = 1.0 / (i_grid_m [i + 1] - i_grid_m [i - 1]);
				}
				oodz [m - 1] = 1.0 / (i_grid_m [m - 1] - i_grid_m [m - 2]);
				
				setup ();
			}
		
			virtual ~background_vertical () {}

			virtual int type () {
				return plan <datatype>::pre;
			}
			
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
				if (*component_flags & z_solve) {
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
				virtual std::shared_ptr <plan <datatype> > _instance (grids::grid <datatype> **grids, datatype **matrices, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
					return std::shared_ptr <plan <datatype> > (new background_vertical <datatype> (*grids [0], *grids [1], alpha, diffusion, explicit_calculate, matrices [0], matrices [1], i_data_in, i_data_out, i_element_flags, i_component_flags));
				}
			};
		};
	} /* diffusion */
} /* plans */

#endif /* end of include guard: DIFFUSION_TWO_D_HPP_YHECX9VS */
