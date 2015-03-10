/*!**********************************************************************
 * \file source_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces/src/two_d
 * 
 * Created by Justin Brown on 2013-12-24.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOURCE_TWO_D_HPP_G9AN5CH6
#define SOURCE_TWO_D_HPP_G9AN5CH6

#include <omp.h>

#include "../explicit_plan.hpp"
#include "../real_plan.hpp"
#include "linalg/exceptions.hpp"
#include "io/parameters.hpp"

namespace plans
{
	template <class datatype>
	class source : public explicit_plan <datatype>
	{
	private:
		using explicit_plan <datatype>::n;
		using explicit_plan <datatype>::ldn;
		using explicit_plan <datatype>::m;
		using explicit_plan <datatype>::data_out;
		
		datatype coeff; //!< The coefficient of the source term
		datatype *data_source; //!< The data pointer for the source data
		
	public:
		/*!**********************************************************************
		 * \copydoc explicit_plan::explicit_plan
		 * 
		 * \param i_coeff The coefficient for the source term
		 * \param i_data_source The data pointer for the source data
		 * 
		 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
		 ************************************************************************/
		source (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype i_coeff, datatype* i_data_source, datatype *i_data_in, datatype *i_data_out, int *i_element_flags, int *i_component_flags) : explicit_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out, i_element_flags, i_component_flags), coeff (i_coeff), data_source (i_data_source) {
			TRACE ("Adding source...");
		}
		
		virtual ~source () {}
		
		/*!**********************************************************************
		 * \copydoc explicit_plan::execute
		 ************************************************************************/
		virtual void execute () {
			TRACE ("Executing source...");
			linalg::matrix_add_scaled (m, ldn, coeff, data_source, data_out, m, m);	
		}
		
		/*!**********************************************************************
		 * \copydoc explicit_plan::factory
		 ************************************************************************/
		class factory : public explicit_plan <datatype>::factory
		{
		private:
			datatype coeff; //!< The coefficient to be used when constructing the plan
			datatype *data_source; //!< The data source to be used when constructing the plan
			
		public:
			/*!**********************************************************************
			 * \param i_coeff The coefficient to be used when constructing the plan
			 * \param i_data_source The data source to be used when constructing the plan
			 ************************************************************************/
			factory (datatype i_coeff, datatype *i_data_source) : coeff (i_coeff), data_source (i_data_source) {}
			
			/*!**********************************************************************
			 * \param i_coeff A YAML::Node scalar object from which the coefficient for the plan should be read
			 * \param i_data_source The data source to be used when constructing the plan
			 * 
			 * If the YAML::Node is not defined, no instance will be created
			 ************************************************************************/
			factory (YAML::Node i_coeff, datatype *i_data_source) : data_source (i_data_source) {
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
					return std::shared_ptr <plans::plan <datatype> > (new source <datatype> (*grids [0], *grids [1], coeff, data_source, i_data_in, i_data_out, i_element_flags, i_component_flags));
				}
				DEBUG ("NOPE");
				return NULL;
			}
		};
	};
} /* plans */

#endif /* end of include guard: SOURCE_TWO_D_HPP_G9AN5CH6 */

