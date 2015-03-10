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
			using real_plan <datatype>::n;
			using real_plan <datatype>::m;
			using real_plan <datatype>::grid_n;
			using real_plan <datatype>::grid_m;
			using real_plan <datatype>::data_in;
			using real_plan <datatype>::data_out;
		
			datatype coeff; //!< The coefficient of the advection term
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
			uniform (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype i_coeff, datatype* i_vel_n, datatype *i_vel_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) : real_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out, i_element_flags, i_component_flags), coeff (-i_coeff), vel_n (i_vel_n), vel_m (i_vel_m), pos_n (&(grid_n [0])), pos_m (&(grid_m [0])) {
				TRACE ("Adding advection...");
				x_vec.resize (n * m);
				x_ptr = &x_vec [0];
				z_vec.resize (n * m);
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
		
			/*!**********************************************************************
			 * \copydoc real_plan::execute
			 ************************************************************************/
			virtual void execute () {
				// Calculate the x component
				linalg::scale (m, 0.0, x_ptr + m * (n - 1));
				linalg::matrix_copy (m, n - 1, data_in + m, x_ptr);
				linalg::matrix_add_scaled (m, n - 1, -1.0, data_in, x_ptr + m);
				linalg::add_scaled (m, 1.0, data_in, x_ptr + m * (n - 1));
				linalg::add_scaled (m, -1.0, data_in + m * (n - 1), x_ptr);
			
				// Calculate the z component
				linalg::scale (n, 0.0, z_ptr + m - 1, m);
				linalg::matrix_copy (m - 1, n, data_in + 1, z_ptr, m, m);
				linalg::matrix_add_scaled (m - 1, n, -1.0, data_in, z_ptr + 1, m, m);
				linalg::add_scaled (n, 1.0, data_in, z_ptr, m, m);
				linalg::add_scaled (n, -1.0, data_in + m - 1, z_ptr + m - 1, m, m);
			
				#pragma omp parallel for
				for (int j = 0; j < m; ++j) {
					for (int i = 0; i < n; ++i) {
						x_ptr [i * m + j] = (vel_n [i * m + j] * x_ptr [i * m + j] * oodx_ptr [i] + vel_m [i * m + j] * z_ptr [i * m + j] * oodz_ptr [j]);
					}
				}
			
				// Scale the whole thing by the coefficient
				linalg::matrix_add_scaled (m, n, coeff, x_ptr, data_out);
			}
		
			/*!**********************************************************************
			 * \copydoc real_plan::factory
			 ************************************************************************/
			class factory : public real_plan <datatype>::factory
			{
			private:
				datatype coeff; //!< The coefficient of the plan to be created
				datatype *vel_n; //!< A pointer to the horizontal component of the velocity
				datatype *vel_m; //!< A pointer to the vertical component of the velocity
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient of the plan to be created
				 * \param i_vel_n A pointer to the horizontal component of the velocity
				 * \param i_vel_m A pointer to the vertical component of the velocity
				 ************************************************************************/
				factory (datatype i_coeff, datatype* i_vel_n, datatype* i_vel_m) : coeff (i_coeff), vel_n (i_vel_n), vel_m (i_vel_m) {}
			
				/*!**********************************************************************
				 * \param i_coeff A YAML::Node scalar to be read into the coefficient of the plan to construct
				 * \param i_vel_n A pointer to the horizontal component of the velocity
				 * \param i_vel_m A pointer to the vertical component of the velocity
				 * 
				 * If the YAML::Node is not defined, no plan will be constructed when instance is called.
				 ************************************************************************/
				factory (YAML::Node i_coeff, datatype* i_vel_n, datatype* i_vel_m) : vel_n (i_vel_n), vel_m (i_vel_m) {
					if (i_coeff.IsDefined ()) {
						coeff = i_coeff.as <datatype> ();
					} else {
						coeff = 0.0;
					}
				}
		
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc real_plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > instance (grids::grid <datatype> **grids, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
					if (coeff) {
						return std::shared_ptr <plans::plan <datatype> > (new uniform <datatype> (*grids [0], *grids [1], coeff, vel_n, vel_m, i_data_in, i_data_out, i_element_flags, i_component_flags));
					}
					return NULL;
				}
			};
		};
	} /* advection */
} /* plans */

#endif /* end of include guard: ADVECTION_TWO_D_HPP_GGR0NN1Q */
