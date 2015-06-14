/*!**********************************************************************
 * \file plans/plans-source/source.hpp
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
		class uniform : public explicit_plan <datatype>
		{
		private:
			using explicit_plan <datatype>::coeff;
			using explicit_plan <datatype>::n;
			using explicit_plan <datatype>::ldn;
			using explicit_plan <datatype>::m;
			using explicit_plan <datatype>::data_out;
		
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
			uniform (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype* i_data_source, datatype *i_data_in, datatype *i_data_out, datatype i_coeff = 1.0, int *i_element_flags = NULL, int *i_component_flags = NULL) : explicit_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out, i_coeff, i_element_flags, i_component_flags), data_source (i_data_source) {
				TRACE ("Adding source...");
				DEBUG (" " << coeff);
			}
		
			virtual ~uniform () {}

			virtual int type () {
				return plan <datatype>::mid;
			}
		
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
				datatype *data_source; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 ************************************************************************/
				factory (datatype *i_data_source, datatype i_coeff = 1.0) : explicit_plan <datatype>::factory (i_coeff), data_source (i_data_source) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc explicit_plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > _instance (grids::grid <datatype> **grids, datatype **matrices, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
					if (coeff) {
						DEBUG (data_source << " " << i_data_in << " " << i_data_out);
						return std::shared_ptr <plans::plan <datatype> > (new uniform <datatype> (*grids [0], *grids [1], data_source, i_data_in, i_data_out, 1.0, i_element_flags, i_component_flags));
					}
					return std::shared_ptr <plans::plan <datatype> > ();
				}
			};
		};
	} /* source */
} /* plans */

#endif /* end of include guard: SOURCE_TWO_D_HPP_G9AN5CH6 */

