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
		template <class datatype>
		class uniform : public real_plan <datatype>
		{
		private:
			using real_plan <datatype>::coeff;
			using real_plan <datatype>::n;
			using real_plan <datatype>::m;
			using real_plan <datatype>::dims;
			using real_plan <datatype>::grid_n;
			using real_plan <datatype>::grid_m;
			using real_plan <datatype>::data_in;
			using real_plan <datatype>::data_out;
		
			datatype *vel_n; //!< A pointer to the horizontal component of the velocity
			datatype *vel_m; //!< A pointer to the vertical component of the velocity
			const datatype *pos_n; //!< A pointer to the horizontal position
			const datatype *pos_m; //!< A pointer to the vertical position
		
			std::vector <datatype> x_vec; //!< A vector to hold the x component of the calculation
			std::vector <datatype> z_vec; //!< A vector to hold the z component of the calculation
			std::vector <datatype> oodx_vec; //!< A vector to hold the 1/dx component of the calculation
			std::vector <datatype> oodz_vec; //!< A vector to hold the 1/dz component of the calculation
		
			datatype *x_ptr; //!< A pointer to x_vec for speed
			datatype *z_ptr; //!< A pointer to z_vec for speed
			datatype *oodx_ptr; //!< A pointer to oodx_vec for speed
			datatype *oodz_ptr; //!< A pointer to oodz_vec for speed
		
		public:
			/*!**********************************************************************
			 * \copydoc real_plan::real_plan
			 * 
			 * \param i_coeff The coefficient of the advection term when it's on the left hand side of the equation
			 * \param i_vel_n A pointer to the horizontal component of the velocity
			 * \param i_vel_m A pointer to the vertical component of the velocity
			 ************************************************************************/
			uniform (grids::variable <datatype> &i_vel_n, grids::variable <datatype> &i_vel_m, grids::variable <datatype> &i_data_in, datatype *i_data_out = NULL, datatype i_coeff = 1.0, int *i_element_flags = NULL, int *i_component_flags = NULL) : real_plan <datatype> (i_data_in, i_data_out, i_coeff, i_element_flags, i_component_flags), vel_n (i_vel_n.ptr ()), vel_m (i_vel_m.ptr ()), pos_n (&(grid_n [0])), pos_m (&(grid_m [0])) {
				TRACE ("Adding advection...");
				x_vec.resize (n * m * dims);
				x_ptr = &x_vec [0];
				z_vec.resize (n * m * dims);
				z_ptr = &z_vec [0];
			
				// Calculate some derivative contributions that will be used frequently
				oodx_vec.resize (n);
				oodx_ptr = &oodx_vec [0];
				oodz_vec.resize (m);
				oodz_ptr = &oodz_vec [0];
			
				// Calculate 1/dx
				oodx_ptr [0] = 0.5 / (pos_n [1] - pos_n [0]);
				for (int i = 1; i < n - 1; ++i) {
					oodx_ptr [i] = 1.0 / (pos_n [i + 1] - pos_n [i - 1]);
				}
				oodx_ptr [n - 1] = 0.5 / (pos_n [n - 1] - pos_n [n - 2]);
			
				// Calculate 1/dz
				oodz_ptr [0] = 1.0 / (pos_m [1] - pos_m [0]);
				for (int i = 1; i < m - 1; ++i) {
					oodz_ptr [i] = 1.0 / (pos_m [i + 1] - pos_m [i - 1]);
				}
				oodz_ptr [m - 1] = 1.0 / (pos_m [m - 1] - pos_m [m - 2]);
			}
	
			virtual ~uniform () {}

			virtual int type () {
				return plan <datatype>::post;
			}
		
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
			}
		
			/*!**********************************************************************
			 * \copydoc real_plan::factory
			 ************************************************************************/
			class factory : public real_plan <datatype>::factory
			{
			private:
				grids::variable <datatype> &vel_n; //!< A pointer to the horizontal component of the velocity
				grids::variable <datatype> &vel_m; //!< A pointer to the vertical component of the velocity
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient of the plan to be created
				 * \param i_vel_n A pointer to the horizontal component of the velocity
				 * \param i_vel_m A pointer to the vertical component of the velocity
				 ************************************************************************/
				factory (grids::variable <datatype> &i_vel_n, grids::variable <datatype> &i_vel_m, datatype i_coeff = -1.0) : real_plan <datatype>::factory (i_coeff), vel_n (i_vel_n), vel_m (i_vel_m) {}
		
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc real_plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
					if (coeff) {
						return std::shared_ptr <plans::plan <datatype> > (new uniform <datatype> (vel_n, vel_m, i_data_in, i_data_out, 1.0, i_element_flags, i_component_flags));
					}
					return std::shared_ptr <plans::plan <datatype> > ();
				}
			};
		};
	} /* advection */
} /* plans */

#endif /* end of include guard: ADVECTION_TWO_D_HPP_GGR0NN1Q */
