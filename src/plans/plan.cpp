/*!***********************************************************************
 * \file plan.hpp
 * Spectral Element
 * 
 * This file provides the abstract base class plan, from which all code 
 * executions should derive.
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "plan.hpp"
#include "source.hpp"

 namespace plans
 {
 	plan::factory_container plan::factory_container::operator+ (grids::variable &var) {
 		return *this + src (var, false);
 	}

 	plan::factory_container plan::factory_container::operator+ (double scalar) {
 		return *this + constant (scalar);
 	}

	plan::factory_container operator+ (std::shared_ptr <plan::factory> i_factory, std::shared_ptr <plan::factory> j_factory) {
		return plan::factory_container (i_factory) + j_factory;
	}

	plan::factory_container operator+ (std::shared_ptr <plan::factory> i_factory, plan::factory_container j_container) {
		return j_container + i_factory;
	}

	plan::factory_container operator- (std::shared_ptr <plan::factory> i_factory) {
		return i_factory * (-1.);
	}

	std::shared_ptr <plan::factory> operator* (std::shared_ptr <plan::factory> i_factory, double scalar) {
		i_factory->coeff *= scalar;
		return i_factory;
	}

	// template <>
	// plan::factory_container operator* (plan::factory_container i_container, double scalar) {
	// 	plan::factory_container container (i_container);
	// 	for (int i = 0; i < (int) i_container.facts.size (); ++i)
	// 	{
	// 		container.facts [i]->coeff *= scalar;
	// 	}
	// 	return container;
	// }

	std::shared_ptr <plan::factory> operator* (std::shared_ptr <plan::factory> i_factory, YAML::Node node) {
		i_factory->coeff *= node.as <double> ();
		return i_factory;
	}

	/**
	 * @copydoc operator*(datatype, plan<datatype>::factory_container)
	 */
	plan::factory_container operator* (double scalar, plan::factory_container i_container) {
		return i_container * scalar;
	}
 }