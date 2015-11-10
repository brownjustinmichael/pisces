/*!**********************************************************************
 * \file plans/source.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOURCE_HPP_206A54EA
#define SOURCE_HPP_206A54EA

#include "plans-source/uniform.hpp"
#include "plans-source/deriv.hpp"
#include "plans-source/heat.hpp"

namespace plans
{
	/**
	 * @brief Shorthand to construct a uniform source term
	 * 
	 * @param data_source A reference to the data to source in the plan
	 * @param dealias Whether to ignore the last third of the data in the horizontal direction
	 * @return A shared pointer the plan factory of the vertical gradient
	 */
	std::shared_ptr <plan::factory> src (grids::variable &data_source, bool dealias = false);

	/**
	 * @brief Shorthand to construct a horizontal gradient term
	 * 
	 * @param data_source A reference to the data over which the derivative should be taken
	 * @return A shared pointer the plan factory of the horizontal gradient
	 */
	std::shared_ptr <plan::factory> grad_x (grids::variable &data_source);

	/**
	 * @brief Shorthand to construct a vertical gradient term
	 * 
	 * @param data_source A reference to the data over which the derivative should be taken
	 * @return A shared pointer the plan factory of the vertical gradient
	 */
	std::shared_ptr <plan::factory> grad_z (grids::variable &data_source);

	/**
	 * @brief Shorthand to construct a constant source term
	 * 
	 * @param coeff The value of the constant source
	 * @return A shared pointer to the plan factory of the constant source
	 */
	std::shared_ptr <plan::factory> constant (double coeff = 1.0);

	/**
	 * @brief Shorthand to construct a viscous heating term
	 * 
	 * @param data_source A multiplier source term in front of the viscous heat term
	 * @param data_x The x component of the velocity
	 * @param data_z The z component of the velocity
	 * @return A shared pointer to the plan factory of the heat term
	 */
	plan::factory_container heat (grids::variable &data_source, grids::variable &data_x, grids::variable &data_z);

	/**
	 * @brief Shorthand to construct a divergence heating term
	 * 
	 * @param data_source A multiplier source term in front of the divergence
	 * @param data_x The x component of the velocity
	 * @param data_z The z component of the velocity
	 * @return A shared pointer to the plan factory of the devergence term
	 */
	plan::factory_container diverge (grids::variable &data_source, grids::variable &data_x, grids::variable &data_z);
} /* plans */

#endif /* end of include guard: SOURCE_HPP_206A54EA */
