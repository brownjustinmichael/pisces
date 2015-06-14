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
	typename explicit_plan <datatype>::factory_container src (datatype *data_source, datatype coeff = 1.0) {
		return typename explicit_plan <datatype>::factory_container (std::shared_ptr <typename explicit_plan <datatype>::factory> (new typename source::uniform <datatype>::factory (data_source, coeff)));
	}
} /* plans */

#endif /* end of include guard: SOURCE_HPP_206A54EA */
