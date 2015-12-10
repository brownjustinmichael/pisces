/*!**********************************************************************
 * \file plans-advection/advection.hpp
 * /Users/justinbrown/Dropbox/pisces/src/two_d
 * 
 * Created by Justin Brown on 2013-12-16.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ADVECTION_TWO_D_HPP_GGR0NN1Q
#define ADVECTION_TWO_D_HPP_GGR0NN1Q

#include "../real_plan.hpp"
#include "io/parameters.hpp"
#include "linalg/utils.hpp"

/*!**********************************************************************
 * \namespace plans::advection
 * 
 * \brief A namespace containing the advection plan extension
 ************************************************************************/
namespace plans
{
	namespace advection
	{
		/*!**********************************************************************
		 * \brief This class calculates the advection component from Cartesian data
		 * 
		 * This term is assumed to be on the left hand side of the equation, so the coefficient's sign is reversed automatically. Please provide the coefficient of the advection term when it appears on the left side of the equation.
		 ************************************************************************/
		class uniform : public real_plan
		{
		private:
			using real_plan::coeff;
			using real_plan::n;
			using real_plan::m;
			using real_plan::dims;
			using real_plan::grid_n;
			using real_plan::grid_m;
			using real_plan::data_in;
			using real_plan::data_out;
		
			double *vel_n; //!< A pointer to the horizontal component of the velocity
			double *vel_m; //!< A pointer to the vertical component of the velocity
			const double *pos_n; //!< A pointer to the horizontal position
			const double *pos_m; //!< A pointer to the vertical position
		
			std::vector <double> x_vec; //!< A vector to hold the x component of the calculation
			std::vector <double> z_vec; //!< A vector to hold the z component of the calculation
		
			double *x_ptr; //!< A pointer to x_vec for speed
			double *z_ptr; //!< A pointer to z_vec for speed
			double *oodx_ptr; //!< A pointer to oodx_vec for speed
			double *oodz_ptr; //!< A pointer to oodz_vec for speed
		
		public:
			/*!**********************************************************************
			 * \copydoc real_plan::real_plan
			 * 
			 * \param i_coeff The coefficient of the advection term when it's on the left hand side of the equation
			 * \param i_vel_n A pointer to the horizontal component of the velocity
			 * \param i_vel_m A pointer to the vertical component of the velocity
			 ************************************************************************/
			uniform (grids::variable &i_vel_n, grids::variable &i_vel_m, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0) : 
			real_plan (i_data_in, i_data_out, i_coeff), vel_n (i_vel_n.ptr ()), vel_m (i_vel_m.ptr ()), pos_n (&(grid_n [0])), pos_m (&(grid_m [0])) {
				TRACE ("Adding advection...");
				x_vec.resize (n * m * dims);
				x_ptr = &x_vec [0];
				z_vec.resize (n * m * dims);
				z_ptr = &z_vec [0];
			
				// Calculate some derivative contributions that will be used frequently
				oodx_ptr = grid_n.get_ood2 ();
				oodz_ptr = grid_m.get_ood2 ();
			}
	
			virtual ~uniform () {}

			/*!**********************************************************************
			 * \copydoc real_plan::execute
			 ************************************************************************/
			virtual void execute () {
				// Calculate the x component
				linalg::scale (m * dims, 0.0, x_ptr + m * (n - 1) * dims);
				linalg::matrix_copy (m * dims, n - 1, data_in + m * dims, x_ptr);
				linalg::matrix_add_scaled (m * dims, n - 1, -1.0, data_in, x_ptr + m * dims);
				linalg::add_scaled (m * dims, 1.0, data_in, x_ptr + m * (n - 1) * dims);
				linalg::add_scaled (m * dims, -1.0, data_in + m * (n - 1) * dims, x_ptr);
			
				// Calculate the z component
				linalg::matrix_scale (dims, n, 0.0, z_ptr + (m - 1) * dims, m * dims);
				linalg::matrix_copy ((m - 1) * dims, n, data_in + dims, z_ptr, m * dims, m * dims);
				linalg::matrix_add_scaled ((m - 1) * dims, n, -1.0, data_in, z_ptr + dims, m * dims, m * dims);
				linalg::matrix_add_scaled (dims, n, 1.0, data_in, z_ptr, m * dims, m * dims);
				linalg::matrix_add_scaled (dims, n, -1.0, data_in + (m - 1) * dims, z_ptr + (m - 1) * dims, m * dims, m * dims);
			
				#pragma omp parallel for
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j < m; ++j) {
						for (int k = 0; k < dims; ++k)
						{
							x_ptr [(i * m + j) * dims + k] = (vel_n [i * m + j] * x_ptr [(i * m + j) * dims + k] * oodx_ptr [i] + vel_m [i * m + j] * z_ptr [i * m + j] * oodz_ptr [j]);
						}
					}
				}
			
				// Scale the whole thing by the coefficient
				linalg::matrix_add_scaled (m * dims, n, coeff, x_ptr, data_out);

				// #pragma omp parallel for
				// for (int i = 1; i < n - 1; ++i)
				// {
				// 	for (int j = 1; j < m - 1; ++j)
				// 	{
				// 		data_out [i * m + j] += coeff * (vel_m [i * m + j + 1] * data_in [i * m + j + 1] - vel_m [i * m + j - 1] * data_in [i * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]);
				// 		data_out [i * m + j] += coeff * (vel_n [(i + 1) * m + j] * data_in [(i + 1) * m + j] - vel_n [(i - 1) * m + j] * data_in [(i - 1) * m + j]) / (pos_n [i + 1] - pos_n [i - 1]);
				// 	}
				// }

				// for (int j = 1; j < m - 1; ++j)
				// {
				// 	data_out [j] += coeff * (vel_m [j + 1] * data_in [j + 1] - vel_m [j - 1] * data_in [j - 1]) / (pos_m [j + 1] - pos_m [j - 1]);
				// 	data_out [j] += coeff * (vel_n [m + j] * data_in [m + j] - vel_n [j] * data_in [j]) / (pos_n [1] - pos_n [0]);
				// 	data_out [(n - 1) * m + j] += coeff * (vel_m [(n - 1) * m + j + 1] * data_in [(n - 1) * m + j + 1] - vel_m [(n - 1) * m + j - 1] * data_in [(n - 1) * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]);
				// 	data_out [(n - 1) * m + j] += coeff * (vel_n [(n - 1) * m + j] * data_in [(n - 1) * m + j] - vel_n [(n - 2) * m + j] * data_in [(n - 2) * m + j]) / (pos_n [n - 1] - pos_n [n - 2]);
				// }
			}

			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public real_plan::factory
			{
			private:
				grids::variable &vel_n; //!< A pointer to the horizontal component of the velocity
				grids::variable &vel_m; //!< A pointer to the vertical component of the velocity
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient of the plan to be created
				 * \param i_vel_n A pointer to the horizontal component of the velocity
				 * \param i_vel_m A pointer to the vertical component of the velocity
				 ************************************************************************/
				factory (grids::variable &i_vel_n, grids::variable &i_vel_m, double i_coeff = 1.0) : real_plan::factory (i_coeff), vel_n (i_vel_n), vel_m (i_vel_m) {}
		
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan > (new uniform (vel_n, vel_m, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan > ();
				}
			};
		};
	} /* advection */
} /* plans */

#endif /* end of include guard: ADVECTION_TWO_D_HPP_GGR0NN1Q */
