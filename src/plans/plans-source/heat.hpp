/*!**********************************************************************
 * \file plans/plans-source/heat.hpp
 * /Users/justinbrown/Dropbox/pisces/src/two_d
 * 
 * Created by Justin Brown on 2013-12-24.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef HEAT_TWO_D_HPP_G9AN5CH6
#define HEAT_TWO_D_HPP_G9AN5CH6

#include <omp.h>

#include "../real_plan.hpp"
#include "linalg/exceptions.hpp"

namespace plans
{
	namespace source
	{
		/*!**********************************************************************
		 * \brief A plan to add viscous heating
		 ************************************************************************/
		template <class datatype>
		class viscous_heat : public real_plan
		{
		private:
			using real_plan::coeff;
			using real_plan::n;
			using real_plan::ldn;
			using real_plan::m;
			using real_plan::dims;
			using real_plan::data_out;
			using real_plan::grid_n;
			using real_plan::grid_m;
		
			datatype *data_source; //!< The data pointer for the source data
			datatype *data_x; //!< The data pointer to the x component of the velocity
			datatype *data_z; //!< The data pointer to the z component of the velocity
			datatype *oodx2; //!< The data pointer to one over the centered x derivative
			datatype *oodz2; //!< The data pointer to one over the centered z derivative
		
		public:
			/**
			 * @copydoc real_plan::real_plan
			 * 
			 * @param i_data_source A reference to the source data to use for the heating
			 * @param i_data_x A reference to the x component of the velocity
			 * @param i_data_z A reference to the z component of the velocity
			 */
			viscous_heat (grids::variable &i_data_source, grids::variable &i_data_x, grids::variable &i_data_z, grids::variable &i_data_in, grids::variable &i_data_out, datatype i_coeff = 1.0) : 
			real_plan (i_data_in, i_data_out, i_coeff), 
			data_source (i_data_source.ptr ()),
			data_x (i_data_x.ptr ()),
			data_z (i_data_z.ptr ()) {
				TRACE ("Adding source...");
				oodx2 = grid_n.get_ood2 ();
				oodz2 = grid_m.get_ood2 ();
			}
		
			virtual ~viscous_heat () {}

			/*!**********************************************************************
			 * \copydoc real_plan::execute
			 ************************************************************************/
			virtual void execute () {
				TRACE ("Executing source...");
				for (int i = 0; i < n; ++i)
				{
					int p1 = (i + 1) % n;
					int m1 = (i - 1 + n) % n;
					int g = 0;
					datatype dudx, dwdx, dudz, dwdz;
					for (int j = 1; j < m - 1; ++j)
					{
						g = i * m + j;
						dudx = (data_x [p1 * m + j] - data_x [m1 * m + j]) * oodx2 [i];
						dudz = (data_x [g + 1] - data_x [g - 1]) * oodz2 [j];
						dwdx = (data_z [p1 * m + j] - data_z [m1 * m + j]) * oodx2 [i];
						dwdz = (data_z [g + 1] - data_z [g - 1]) * oodz2 [j];

						data_out [g] += coeff * data_source [g] * dudx * (4. / 3. * dudx - 2. / 3. * dwdz);
						data_out [g] += coeff * data_source [g] * dudz * (dudz + dwdx);

						data_out [g] += coeff * data_source [g] * dwdz * (4. / 3. * dwdz - 2. / 3. * dudx);
						data_out [g] += coeff * data_source [g] * dwdx * (dwdx + dudz);
					}
				}
			}
		
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public real_plan::factory
			{
			private:
				grids::variable &data_source; //!< The data source to be used when constructing the plan
				grids::variable &data_x; //!< The data source to be used when constructing the plan
				grids::variable &data_z; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 * \param i_data_x The x component of the velocity data to be used when constructing the plan
				 * \param i_data_z The z component of the velocity data source to be used when constructing the plan
				 ************************************************************************/
				factory (grids::variable &i_data_source, grids::variable &i_data_x, grids::variable &i_data_z, datatype i_coeff = 1.0) : 
				real_plan::factory (i_coeff), 
				data_source (i_data_source),
				data_x (i_data_x),
				data_z (i_data_z) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (datatype **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan > (new viscous_heat <datatype> (data_source, data_x, data_z, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan > ();
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to account for temperature change due to expansion and contraction
		 ************************************************************************/
		template <class datatype>
		class divergence : public real_plan
		{
		private:
			using real_plan::coeff;
			using real_plan::n;
			using real_plan::ldn;
			using real_plan::m;
			using real_plan::dims;
			using real_plan::data_out;
			using real_plan::grid_n;
			using real_plan::grid_m;
		
			datatype *data_source; //!< The data pointer for the source data
			datatype *data_x; //!< The data pointer to the x component of the velocity
			datatype *data_z; //!< The data pointer to the z component of the velocity
			datatype *oodx2; //!< The data pointer to one over the centered x derivative
			datatype *oodz2; //!< The data pointer to one over the centered z derivative
		
		public:
			/**
			 * @copydoc real_plan::real_plan
			 * 
			 * @param i_data_source A reference to the source data to use for the heating
			 * @param i_data_x A reference to the x component of the velocity
			 * @param i_data_z A reference to the z component of the velocity
			 */
			divergence (grids::variable &i_data_source, grids::variable &i_data_x, grids::variable &i_data_z, grids::variable &i_data_in, grids::variable &i_data_out, datatype i_coeff = 1.0) : 
			real_plan (i_data_in, i_data_out, i_coeff), 
			data_source (i_data_source.ptr ()),
			data_x (i_data_x.ptr ()),
			data_z (i_data_z.ptr ()) {
				TRACE ("Adding source...");
				oodx2 = grid_n.get_ood2 ();
				oodz2 = grid_m.get_ood2 ();
			}
		
			virtual ~divergence () {}

			/*!**********************************************************************
			 * \copydoc real_plan::execute
			 ************************************************************************/
			virtual void execute () {
				TRACE ("Executing source...");
				for (int i = 0; i < n; ++i)
				{
					int p1 = (i + 1) % n;
					int m1 = (i - 1 + n) % n;
					int g = 0;
					for (int j = 1; j < m - 1; ++j)
					{
						g = i * m + j;
						data_out [g] += coeff * data_source [g] * ((data_x [p1 * m + j] - data_x [m1 * m + j]) * oodx2 [i] + (data_z [g + 1] - data_z [g - 1]) * oodz2 [j]);
					}
				}
			}
		
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public real_plan::factory
			{
			private:
				grids::variable &data_source; //!< The data source to be used when constructing the plan
				grids::variable &data_x; //!< The data source to be used when constructing the plan
				grids::variable &data_z; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 * \param i_data_x The x component of the velocity data to be used when constructing the plan
				 * \param i_data_z The z component of the velocity data source to be used when constructing the plan
				 ************************************************************************/
				factory (grids::variable &i_data_source, grids::variable &i_data_x, grids::variable &i_data_z, datatype i_coeff = 1.0) : 
				real_plan::factory (i_coeff), 
				data_source (i_data_source),
				data_x (i_data_x),
				data_z (i_data_z) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (datatype **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan > (new divergence <datatype> (data_source, data_x, data_z, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan > ();
				}
			};
		};
	} /* source */
} /* plans */

#endif /* end of include guard: SOURCE_TWO_D_HPP_G9AN5CH6 */

