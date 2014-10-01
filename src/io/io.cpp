/*!***********************************************************************
 * \file io.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "io.hpp"

namespace io
{	
	YAML::Node parameters::operator [] (std::string key) {
		std::istringstream ss (key);
		std::string token;
		std::getline (ss, token, '.');
		std::vector <YAML::Node> nodes;
		nodes.push_back (YAML::Node::operator [] (token));
		while (std::getline (ss, token, '.')) {
			nodes.push_back (nodes [(int) nodes.size () - 1] [token]);
		}
		return nodes [(int) nodes.size () - 1];
	}
}	

