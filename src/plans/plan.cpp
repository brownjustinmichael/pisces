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

 namespace plans
 {
	typename plan <double>::factory_container operator+ (std::shared_ptr <typename plan <double>::factory> i_factory, std::shared_ptr <typename plan <double>::factory> j_factory) {
		return plan <double>::factory_container (i_factory) + j_factory;
	}

	typename plan <double>::factory_container operator- (std::shared_ptr <typename plan <double>::factory> i_factory, std::shared_ptr <typename plan <double>::factory> j_factory) {
		return plan <double>::factory_container (i_factory) - j_factory;
	}

	template <>
	typename plan <double>::factory_container operator* (typename plan <double>::factory_container i_container, double scalar) {
		typename plan <double>::factory_container container (i_container);
		for (int i = 0; i < (int) i_container.facts.size (); ++i)
		{
			container.facts [i]->coeff *= scalar;
		}
		return container;
	}

	typename std::shared_ptr <typename plan <double>::factory> operator* (std::shared_ptr <typename plan <double>::factory> i_factory, YAML::Node node) {
		i_factory->coeff *= node.as <double> ();
		return i_factory;
	}

	template <>
	typename plan <double>::factory_container operator* (double scalar, typename plan <double>::factory_container i_container) {
		return i_container * scalar;
	}
 }