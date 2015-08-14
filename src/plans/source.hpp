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
#include "plans-source/heat.hpp"

namespace plans
{
	template <class datatype>
	std::shared_ptr <typename plan <datatype>::factory> src (grids::variable <datatype> &data_source, datatype coeff = 1.0) {
		return std::shared_ptr <typename explicit_plan <datatype>::factory> (new typename source::uniform <datatype>::factory (data_source, coeff));
	}

	template <class datatype>
	std::shared_ptr <typename plan <datatype>::factory> grad_x (grids::variable <datatype> &data_source, datatype coeff = 1.0) {
		return std::shared_ptr <typename explicit_plan <datatype>::factory> (new typename source::uniform_grad_x <datatype>::factory (data_source, coeff));
	}

	template <class datatype>
	std::shared_ptr <typename plan <datatype>::factory> grad_z (grids::variable <datatype> &data_source, datatype coeff = 1.0) {
		return std::shared_ptr <typename explicit_plan <datatype>::factory> (new typename source::uniform_grad_z <datatype>::factory (data_source, coeff));
	}

	template <class datatype>
	std::shared_ptr <typename plan <datatype>::factory> z_src (datatype *data_source, datatype coeff = 1.0) {
		return std::shared_ptr <typename real_plan <datatype>::factory> (new typename source::z_src <datatype>::factory (data_source, coeff));
	}

	template <class datatype>
	std::shared_ptr <typename plan <datatype>::factory> constant (datatype coeff = 1.0) {
		return std::shared_ptr <typename explicit_plan <datatype>::factory> (new typename source::constant <datatype>::factory (coeff));
	}

	template <class datatype>
	typename plan <datatype>::factory_container heat (grids::variable <datatype> &data_source, grids::variable <datatype> &data_x, grids::variable <datatype> &data_z, datatype coeff = 1.0) {
		return typename plan <datatype>::factory_container (std::shared_ptr <typename real_plan <datatype>::factory> (new typename source::viscous_heat <datatype>::factory (data_source, data_x, data_z, coeff)));
	}

	template <class datatype>
	typename plan <datatype>::factory_container diverge (grids::variable <datatype> &data_source, grids::variable <datatype> &data_x, grids::variable <datatype> &data_z, datatype coeff = 1.0) {
		return typename plan <datatype>::factory_container (std::shared_ptr <typename real_plan <datatype>::factory> (new typename source::divergence <datatype>::factory (data_source, data_x, data_z, coeff)));
	}

	template <class datatype>
	std::shared_ptr <typename plan <datatype>::factory> pressure_grad_1d (grids::variable <datatype> &i_data_top, grids::variable <datatype> &i_data_bot, datatype *i_data_grad, datatype coeff = 1.0) {
		return std::shared_ptr <typename explicit_plan <datatype>::factory> (new typename source::pressure_grad_1d <datatype>::factory (i_data_top, i_data_bot, i_data_grad, coeff));
	}

} /* plans */

#endif /* end of include guard: SOURCE_HPP_206A54EA */
