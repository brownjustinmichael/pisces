/*!**********************************************************************
 * \file plans/source.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOURCE_HPP_206A54EA
#define SOURCE_HPP_206A54EA

#include "plans-source/source.hpp"

namespace plans
{
	template <class datatype>
	typename explicit_plan <datatype>::factory_container src (grids::variable <datatype> &data_source, datatype coeff = 1.0) {
		return typename explicit_plan <datatype>::factory_container (std::shared_ptr <typename explicit_plan <datatype>::factory> (new typename source::uniform <datatype>::factory (data_source, coeff)));
	}

	template <class datatype>
	typename explicit_plan <datatype>::factory_container constant (datatype coeff = 1.0) {
		return typename explicit_plan <datatype>::factory_container (std::shared_ptr <typename explicit_plan <datatype>::factory> (new typename source::constant <datatype>::factory (coeff)));
	}

	template <class datatype>
	typename explicit_plan <datatype>::factory_container pressure_grad_1d (grids::variable <datatype> &i_data_top, grids::variable <datatype> &i_data_bot, datatype *i_data_grad, datatype coeff = 1.0) {
		return typename explicit_plan <datatype>::factory_container (std::shared_ptr <typename explicit_plan <datatype>::factory> (new typename source::pressure_grad_1d <datatype>::factory (i_data_top, i_data_bot, i_data_grad, coeff)));
	}

} /* plans */

#endif /* end of include guard: SOURCE_HPP_206A54EA */
