/*!**********************************************************************
 * \file parameters.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2015-01-27.
 * Copyright 2015 Justin Brown. All rights reserved.
 ************************************************************************/

#include "parameters.hpp"
#include <iostream>

namespace io
{
	YAML::Node parameters::copy (const YAML::Node &from, const YAML::Node to) {
		YAML::Node copyNode;
		if (from.IsScalar ()) {
			return from;
		}
		if (to.IsScalar ()) {
			return to;
		}
		if (to.IsMap ()) {
			for (YAML::const_iterator iter = to.begin (); iter != to.end (); ++iter) {
				copyNode [(iter->first).as <std::string> ()] = iter->second;
			}
		}
		if (from.IsMap ()) {
			for (YAML::const_iterator iter = from.begin (); iter != from.end (); ++iter) {
				if (copyNode [(iter->first).as <std::string> ()].IsDefined ()) {
					copyNode [(iter->first).as <std::string> ()] = copy (iter->second, copyNode [(iter->first).as <std::string> ()]);
				} else {
					copyNode [(iter->first).as <std::string> ()] = iter->second;
				}
			}
		}
		return copyNode;
	}
	
	YAML::Node parameters::operator[] (std::string key) {
		std::istringstream ss (key);
		std::string token;
		std::getline (ss, token, '.');
		std::vector <std::string> tokens;
		std::vector <YAML::Node> nodes;
		tokens.push_back (token);
		nodes.push_back (YAML::Node::operator [] (token));
		while (std::getline (ss, token, '.')) {
			tokens.push_back (token);
			nodes.push_back (nodes [(int) nodes.size () - 1] [token]);
		}
		int i = (int) nodes.size () - 1;
		std::string inner_tokens = "";
		if (!nodes [i].IsDefined ()) {
			YAML::Node current;
			for (int j = i - 1; j >= 0; --j) {
				inner_tokens = tokens [j + 1] + inner_tokens;
				current = nodes [j];
				if (current.IsDefined () && current ["link"].IsDefined ()) {
					if (parameters::operator[] (current ["link"].as <std::string> () + "." + inner_tokens).IsDefined ()) {
						return parameters::operator[] (current ["link"].as <std::string> () + "." + inner_tokens);
					}
				}
				inner_tokens = "." + inner_tokens;
			}
		}
		return nodes [i];
	}
} /* io */
