/*!**********************************************************************
 * \file source_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces/src/two_d
 * 
 * Created by Justin Brown on 2013-12-24.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOURCE_TWO_D_HPP_G9AN5CH6
#define SOURCE_TWO_D_HPP_G9AN5CH6

#include "../explicit_plan.hpp"
#include "../real_plan.hpp"
#include "linalg/exceptions.hpp"
#include <omp.h>
#include "io/parameters.hpp"

namespace plans
{
	template <class datatype>
	class source : public explicit_plan <datatype>
	{
	public:
		source (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, datatype i_coeff, datatype* i_data_source, datatype *i_data_in, datatype *i_data_out, int *i_element_flags, int *i_component_flags) :
		explicit_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out, i_element_flags, i_component_flags),
		coeff (i_coeff),
		data_source (i_data_source) {
			TRACE ("Adding source...");
		}
	
		virtual ~source () {}
	
		virtual void execute () {
			TRACE ("Executing source...");
			linalg::matrix_add_scaled (m, ldn, coeff, data_source, data_out, m, m);	
		}
	
		class factory : public explicit_plan <datatype>::factory
		{
		private:
			datatype coeff;
			datatype *data_source;
		public:
			factory (datatype i_coeff, datatype *i_data_source) : coeff (i_coeff), data_source (i_data_source) {}
			
			factory (YAML::Node i_coeff, datatype *i_data_source) : data_source (i_data_source) {
				if (i_coeff.IsDefined ()) {
					coeff = i_coeff.as <datatype> ();
				} else {
					coeff = 0.0;
				}
			}
			
			virtual ~factory () {}

			virtual std::shared_ptr <plans::plan <datatype> > instance (plans::grid <datatype> **grids, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
				if (coeff) {
					return std::shared_ptr <plans::plan <datatype> > (new source <datatype> (*grids [0], *grids [1], coeff, data_source, i_data_in, i_data_out, i_element_flags, i_component_flags));
				}
				DEBUG ("NOPE");
				return NULL;
			}
		};

	private:
		using explicit_plan <datatype>::n;
		using explicit_plan <datatype>::ldn;
		using explicit_plan <datatype>::m;
		using explicit_plan <datatype>::data_out;
		datatype coeff;
		datatype *data_source;
	};
} /* plans */

#endif /* end of include guard: SOURCE_TWO_D_HPP_G9AN5CH6 */

