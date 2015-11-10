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
	/**
	 * @brief Shorthand to include both directions of diffusion plans
	 * 
	 * @param alpha The implicit parameter (1.0 is fully implicit, 0.0 is fully explicit)
	 * 
	 * @return A factory container containing the diffusion plans
	 */
	template <class datatype>
	typename plan::factory_container diff (datatype alpha = 1.0) {
		return typename implicit_plan::factory_container (std::shared_ptr <typename implicit_plan::factory> (new typename diffusion::vertical <datatype>::factory (1.0, alpha))) + 
		typename implicit_plan::factory_container (std::shared_ptr <typename implicit_plan::factory> (new typename diffusion::horizontal <datatype>::factory (1.0, alpha)));
	}

	/**
	 * @brief Shorthand to include both directions of background diffusion plans
	 * 
	 * @param i_diffusion A pointer to the background diffusion array (one-dimensional in the z-direction)
	 * @param alpha The implicit parameter (1.0 is fully implicit, 0.0 is fully explicit)
	 * 
	 * @return A factory container containing the diffusion plans
	 */
	template <class datatype>
	typename plan::factory_container bg_diff (datatype *i_diffusion, datatype alpha = 1.0) {
		return typename implicit_plan::factory_container (std::shared_ptr <typename plan::factory> (new typename diffusion::background_vertical <datatype>::factory (alpha, i_diffusion))) + 
		typename implicit_plan::factory_container (std::shared_ptr <typename plan::factory> (new typename diffusion::background_horizontal <datatype>::factory (alpha, i_diffusion)));
	}

	/**
	 * @brief Shorthand to include density-weighted diffusion plans
	 * 
	 * @param density A referencec to the density variable
	 * @param alpha The implicit parameter (1.0 is fully implicit, 0.0 is fully explicit)
	 * 
	 * @return A factory container containing the diffusion plans
	 */
	template <class datatype>
	typename plan::factory_container density_diff (grids::variable &density, datatype alpha = 1.0) {
		return diff <datatype> (alpha) +
		typename plan::factory_container (std::shared_ptr <typename plan::factory> (new typename diffusion::variable_diffusion <datatype>::factory (density, 1.0)));
	}

	/**
	 * @brief Shorthand to include stress terms for the horizontal component of the velocity
	 * 
	 * @param density A reference to the density variable
	 * @param data_other A reference to the other component of the velocity
	 * @return A shared pointer to the factory produced by this
	 */
	std::shared_ptr <typename plan::factory> horizontal_stress (grids::variable &density, grids::variable &data_other);

	/**
	 * @brief Shorthand to include stress terms for the vertical component of the velocity
	 * 
	 * @param density A reference to the density variable
	 * @param data_other A reference to the other component of the velocity
	 * @return A shared pointer to the factory produced by this
	 */
	std::shared_ptr <typename plan::factory> vertical_stress (grids::variable &density, grids::variable &data_other);
} /* plans */

#endif /* end of include guard: DIFFUSION_HPP_C1935DDC */
