/*!**********************************************************************
 * \file plans/plans-source/local_source.hpp
 * /Users/justinbrown/Dropbox/pisces/src/two_d
 * 
 * Created by Justin Brown on 2013-12-24.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef LOCAL_SOURCE_TWO_D_HPP_G9AN5CH6
#define LOCAL_SOURCE_TWO_D_HPP_G9AN5CH6

#include <omp.h>

#include "../explicit_plan.hpp"
#include "../real_plan.hpp"
#include "linalg/exceptions.hpp"
#include "io/parameters.hpp"

/*!**********************************************************************
 * \namespace plans::source
 * 
 * \brief A namespace containing the source plan extension
 ************************************************************************/
namespace plans
{
	namespace source
	{
		/*!**********************************************************************
		 * \brief A plan to add a source term to an equation
		 ************************************************************************/
		template <class datatype>
		class local : public explicit_plan <datatype>
		{
		private:
			using explicit_plan <datatype>::n;
			using explicit_plan <datatype>::ldn;
			using explicit_plan <datatype>::m;
			using explicit_plan <datatype>::data_out;
		
			datatype coeff; //!< The coefficient of the source term
			datatype z_begin;
			datatype z_end;
			int j_begin, j_end;
			datatype *data_source; //!< The data pointer for the source data
		
		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * \param i_data_source The data pointer for the source data
			 * 
			 * In this plan, data_source is not used in lieu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			local (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype i_coeff, datatype i_z_begin, datatype i_z_end, datatype* i_data_source, datatype *i_data_in, datatype *i_data_out, int *i_element_flags, int *i_component_flags) : explicit_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out, i_element_flags, i_component_flags), coeff (i_coeff), data_source (i_data_source) {
				TRACE ("Adding source...");

				setup ();
			}
		
			virtual ~local () {}
			
			virtual void setup () {
				for (int j = 0; j < m; ++j)
				{
					if (grid_m [j] > z_begin)
					{
						j_begin = j;
						break;
					}
				}
				for (int j = m - 1; j <= 0; --j)
				{
					if (grid_m [j] < z_end)
					{
						j_end = j;
						break;
					}
				}
			}

			/*!**********************************************************************
			 * \copydoc explicit_plan::execute
			 ************************************************************************/
			virtual void execute () {
				TRACE ("Executing source...");
				linalg::matrix_add_scaled (j_end - j_begin, ldn, coeff, data_source + j_begin, data_out + j_begin, m, m);	
			}
		
			/*!**********************************************************************
			 * \copydoc explicit_plan::factory
			 ************************************************************************/
			class factory : public explicit_plan <datatype>::factory
			{
			private:
				datatype coeff; //!< The coefficient to be used when constructing the plane
				datatype z_begin, z_end;
				datatype *data_source; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 ************************************************************************/
				factory (datatype i_coeff, datatype i_z_begin, datatype i_z_end, datatype *i_data_source) : coeff (i_coeff), z_begin (i_z_begin), z_end (i_z_end), data_source (i_data_source) {}
			
				/*!**********************************************************************
				 * \param i_coeff A YAML::Node scalar object from which the coefficient for the plan should be read
				 * \param i_data_source The data source to be used when constructing the plan
				 * 
				 * If the YAML::Node is not defined, no instance will be created
				 ************************************************************************/
				factory (YAML::Node i_coeff, datatype i_z_begin, datatype i_z_end, datatype *i_data_source) : z_begin (i_z_begin), z_end (i_z_end), data_source (i_data_source) {
					if (i_coeff.IsDefined ()) {
						coeff = i_coeff.as <datatype> ();
					} else {
						coeff = 0.0;
					}
				}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc explicit_plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > instance (grids::grid <datatype> **grids, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
					if (coeff) {
						return std::shared_ptr <plans::plan <datatype> > (new uniform <datatype> (*grids [0], *grids [1], coeff, z_begin, z_end, data_source, i_data_in, i_data_out, i_element_flags, i_component_flags));
					}
					return std::shared_ptr <plans::plan <datatype> > ();
				}
			};
		};
	} /* source */
} /* plans */

#endif /* end of include guard: LOCAL_SOURCE_TWO_D_HPP_G9AN5CH6 */

