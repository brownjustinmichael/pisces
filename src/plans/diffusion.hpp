/*!**********************************************************************
 * \file diffusion.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DIFFUSION_HPP_C1935DDC
#define DIFFUSION_HPP_C1935DDC

#include "plans-diffusion/horizontal_diffusion.hpp"
#include "plans-diffusion/vertical_diffusion.hpp"
#include "plans-diffusion/variable_diffusion.hpp"

#include "implicit_plan.hpp"
#include "io/parameters.hpp"

namespace plans
{
	template <class datatype>
	typename plan <datatype>::factory_container diff (datatype alpha = 1.0, datatype coeff = 1.0) {
		return typename implicit_plan <datatype>::factory_container (std::shared_ptr <typename implicit_plan <datatype>::factory> (new typename diffusion::vertical <datatype>::factory (coeff, alpha))) + 
		typename implicit_plan <datatype>::factory_container (std::shared_ptr <typename implicit_plan <datatype>::factory> (new typename diffusion::horizontal <datatype>::factory (coeff, alpha)));
	}

	template <class datatype>
	typename plan <datatype>::factory_container density_diff (grids::variable <datatype> &density, datatype alpha = 1.0, datatype coeff = 1.0) {
		return diff <datatype> (alpha, coeff) +
		typename plan <datatype>::factory_container (std::shared_ptr <typename plan <datatype>::factory> (new typename diffusion::variable_diffusion <datatype>::factory (density, coeff)));
	}

	template <class datatype>
	std::shared_ptr <typename plan <datatype>::factory> horizontal_stress (grids::variable <datatype> &density, grids::variable <datatype> &data_other, datatype coeff = 1.0) {
		return std::shared_ptr <typename explicit_plan <datatype>::factory> (new typename diffusion::horizontal_stress <datatype>::factory (density, data_other, coeff));
	}

	template <class datatype>
	std::shared_ptr <typename plan <datatype>::factory> vertical_stress (grids::variable <datatype> &density, grids::variable <datatype> &data_other, datatype coeff = 1.0) {
		return std::shared_ptr <typename explicit_plan <datatype>::factory> (new typename diffusion::vertical_stress <datatype>::factory (density, data_other, coeff));
	}
} /* plans */

#endif /* end of include guard: DIFFUSION_HPP_C1935DDC */
